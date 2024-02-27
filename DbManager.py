import os, math

if __name__ == "__main__":
  print("Standalone operation not available")
  exit()

C_DB_HEADER_SIZE_IN_BYTES    : int = 28
C_TABLE_HEADER_SIZE_IN_BYTES : int =  8
C_APT_ROWSIZE_IN_BYTES       : int = 32
C_VHF_ROWSIZE_IN_BYTES       : int = 40
C_NDB_ROWSIZE_IN_BYTES       : int = 28
C_WPT_ROWSIZE_IN_BYTES       : int = 16

C_EXC_NO_ERROR        : int = -0
C_EXC_INVALID_UNICODE : int = -1
C_EXC_INDEX_ERROR     : int = -2



def SemiCirc2Degrees(SemiC : int) -> float:
  output : float = SemiC * (180.0/math.pow(2,31))
  return output

def LongWord2Int(Chunk : bytes) -> int:
  MyInt = int.from_bytes(Chunk,byteorder="little")
  if MyInt >= math.pow(2,31):
    return MyInt - math.pow(2,32) + 1
  else:
    return MyInt
  
class Record:
  ID     : int
  ICAO   : str
  Lat    : float
  Lon    : float
  Region : str

  def __init__(self) -> None:
    pass

class Database:
  CRC     : int
  Name    : str
  TabSize : dict = {str:int}
  Storage : list[Record]

  def __init__(self) -> None:
    self.Storage = list()
    self.CRC = 0

  def Search(self, key: str) -> list[Record]:
    output : list[Record] = []
    for element in self.Storage:
      if element.ICAO == key:
        output.append(element)
    return output
  
  def Load(self, DbPath : str) -> int:
    Index : int = 0
    ErrorCode : int = C_EXC_NO_ERROR
    with open(DbPath, mode = 'rb') as DbFile:
      Chunk = DbFile.read(C_DB_HEADER_SIZE_IN_BYTES) #Db Header
      ErrorCode = self.ParseDbHeader(Chunk)
      if ErrorCode != 0:
        return ErrorCode
      for index in range(4):
        Chunk = DbFile.read(C_TABLE_HEADER_SIZE_IN_BYTES) #A table header
        TableType = self.ParseTableHeader(Chunk)
        if TableType == 1:
          for index in range(0,self.TabSize["APT"]):
            Chunk = DbFile.read(C_APT_ROWSIZE_IN_BYTES)
            self.Storage.append(ReadAptRecord(AptRow=Chunk))
            if type(self.Storage[-1]) != Record:
              return C_EXC_INVALID_UNICODE
        elif TableType == 2:
          for index in range(0,self.TabSize["VHF"]):
            Chunk = DbFile.read(C_VHF_ROWSIZE_IN_BYTES)
            self.Storage.append(ReadVhfRecord(VhfRow=Chunk))
            if type(self.Storage[-1]) != Record:
              return C_EXC_INVALID_UNICODE
        elif TableType == 3:
          for index in range(0,self.TabSize["NDB"]):
            Chunk = DbFile.read(C_NDB_ROWSIZE_IN_BYTES)
            self.Storage.append(ReadNdbRecord(NdbRow=Chunk))
            if type(self.Storage[-1]) != Record:
              return C_EXC_INVALID_UNICODE
        elif TableType == 4:
          for index in range(0,self.TabSize["WPT"]):
            Chunk = DbFile.read(C_WPT_ROWSIZE_IN_BYTES)
            #self.Storage.append(ReadAptRecord(NdbRow=Chunk))
        else:
          return TableType
    return C_EXC_NO_ERROR

  def ParseDbHeader(self, HeaderChunk : bytes) -> int:
    try:
      SubChunk = HeaderChunk[0:12]
      self.Name = SubChunk.decode(encoding="ascii")
      SubChunk = HeaderChunk[12:15]
      self.CRC = int.from_bytes(SubChunk,byteorder="little")
    except IndexError:
      return C_EXC_INDEX_ERROR
    except UnicodeDecodeError:
      return C_EXC_INVALID_UNICODE
    return C_EXC_NO_ERROR
  
  def ParseTableHeader(self, HeaderChunk: bytes) -> int:
    try:
      TableType : int = 0
      SubChunk = HeaderChunk[0:3]
      TableType   = int.from_bytes(SubChunk,byteorder="little")
      SubChunk = HeaderChunk[4:8]
      TabSize = int.from_bytes(SubChunk,byteorder="little")
    except IndexError:
      return C_EXC_INDEX_ERROR
    except UnicodeDecodeError:
      return C_EXC_INVALID_UNICODE
    if TableType == 1:
      self.TabSize["APT"] = TabSize
    elif TableType == 2:
      self.TabSize["VHF"] = TabSize
    elif TableType == 3:
      self.TabSize["NDB"] = TabSize
    elif TableType == 4:
      self.TabSize["WPT"] = TabSize
    else:
      return C_EXC_INDEX_ERROR
    return TableType
  
def ReadAptRecord(AptRow : bytes) -> Record:
  Aptrecord = Record()
  Aptrecord.ICAO = AptRow[0:8].decode(encoding="ascii")
  Aptrecord.ICAO = Aptrecord.ICAO[0:Aptrecord.ICAO.index("\0")]
  Aptrecord.Region = AptRow[8:12].decode(encoding="ascii")
  Aptrecord.Region = Aptrecord.Region[0:Aptrecord.Region.index("\0")]
  Aptrecord.Lat  = SemiCirc2Degrees(LongWord2Int(Chunk=AptRow[20:24]))
  Aptrecord.Lon  = SemiCirc2Degrees(LongWord2Int(Chunk=AptRow[24:28]))
  return Aptrecord

def ReadVhfRecord(VhfRow : bytes) -> Record:
  Vhfrecord = Record()
  Vhfrecord.ICAO = VhfRow[0:8].decode(encoding="ascii")
  Vhfrecord.ICAO = Vhfrecord.ICAO[0:Vhfrecord.ICAO.index("\0")]
  Vhfrecord.Region = VhfRow[8:12].decode(encoding="ascii")
  Vhfrecord.Region = Vhfrecord.Region[0:Vhfrecord.Region.index("\0")]
  Vhfrecord.Lat  = SemiCirc2Degrees(LongWord2Int(Chunk=VhfRow[16:20]))
  Vhfrecord.Lon  = SemiCirc2Degrees(LongWord2Int(Chunk=VhfRow[20:24]))
  return Vhfrecord

def ReadNdbRecord(NdbRow : bytes) -> Record | int:
  try:
    Ndbrecord = Record()
    Ndbrecord.ICAO = NdbRow[0:8].decode(encoding="ascii")
    Ndbrecord.ICAO = Ndbrecord.ICAO[0:Ndbrecord.ICAO.index("\0")]
    Ndbrecord.Region = NdbRow[8:12].decode(encoding="ascii")
    Ndbrecord.Region = Ndbrecord.Region[0:Ndbrecord.Region.index("\0")]
    Ndbrecord.Lat  = SemiCirc2Degrees(LongWord2Int(Chunk=NdbRow[12:16]))
    Ndbrecord.Lon  = SemiCirc2Degrees(LongWord2Int(Chunk=NdbRow[16:20]))
  except UnicodeDecodeError as BadNdbRecord:
    return 1
  return Ndbrecord