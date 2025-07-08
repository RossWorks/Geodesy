#!/usr/bin/env python3
import tkinter
import Geodesy
import math
import DbManager
from tkinter import messagebox, filedialog, ttk

padValue = 10
CentralWidth = 10
ComputationMethodsList = ["Vincenty", "Sodano", "Sphere", "Flat Earth"]

DefaultFontTuple = ("B612 mono", 15, "normal")
NumericFontTuple = ("B612", 15, "normal")

MyDatabase = DbManager.Database()

SearchResult = list()
SearchResult2 = list()


def LoadDb():
    ErrorCode: int = DbManager.C_EXC_NO_ERROR
    DbFileName: str
    DbFileName = filedialog.askopenfilename()
    if not isinstance(DbFileName, str):
        return
    ErrorCode = MyDatabase.Load(DbPath=DbFileName)
    if ErrorCode != DbManager.C_EXC_NO_ERROR:
        messagebox.showerror(title="Error", message="Invalid DB")
        CmdSearch1.config(state="disabled")
        CmdSearch2.config(state="disabled")
        return
    VhfCount.set(value=MyDatabase.TabSize["VHF"])
    AptCount.set(value=MyDatabase.TabSize["APT"])
    NdbCount.set(value=MyDatabase.TabSize["NDB"])
    LblDbVhfCount.config(text="VHF:" + str(VhfCount.get()).ljust(5, " "))
    LblDbAptCount.config(text="APT:" + str(AptCount.get()).ljust(5, " "))
    LblDbNdbCount.config(text="NDB:" + str(NdbCount.get()).ljust(5, " "))
    FrameDbInfo.config(text=MyDatabase.Name)
    CmdSearch1.config(state="normal")
    CmdSearch2.config(state="normal")
    home.update()


def ConfirmSelection():
    SelectedIndex: int = 0
    SelectedIndex = SearchResultBox.curselection()[0]
    TxtStartLat.delete(0, tkinter.END)
    TxtStartLat.insert(0, str(SearchResult[SelectedIndex].Lat))
    TxtStartLon.delete(0, tkinter.END)
    TxtStartLon.insert(0, str(SearchResult[SelectedIndex].Lon))
    MultipleSel.destroy()


def ConfirmSelection2():
    SelectedIndex: int = 0
    SelectedIndex = SearchResultBox.curselection()[0]
    TxtDestLat.delete(0, tkinter.END)
    TxtDestLat.insert(0, str(SearchResult2[SelectedIndex].Lat))
    TxtDestLon.delete(0, tkinter.END)
    TxtDestLon.insert(0, str(SearchResult2[SelectedIndex].Lon))
    MultipleSel.destroy()


def Inverse_Imp():
    Origin = Geodesy.Point()
    Dest = Geodesy.Point()
    Route = Geodesy.Route()
    TxtStartDist.delete(0, tkinter.END)
    TxtStartBear.delete(0, tkinter.END)
    TxtFinalAz.delete(0, tkinter.END)
    try:
        Origin.Latitude = math.radians(float(TxtStartLat.get()))
        Origin.Longitude = math.radians(float(TxtStartLon.get()))
        Dest.Latitude = math.radians(float(TxtDestLat.get()))
        Dest.Longitude = math.radians(float(TxtDestLon.get()))
    except ValueError:
        TxtStartBear.delete(0, tkinter.END)
        TxtStartBear.insert(0, "INVALID PARAM")
        messagebox.showerror(
            title="Inverse distance",
            message="Check origin and destination point information",
        )
        return
    Method = ListBoxMethod.get()
    if Method == ComputationMethodsList[0]:  # Vincenty
        Route = Geodesy.InverseVincenty(
            OriginPoint=Origin, DestinationPoint=Dest, tol=1e-24
        )
    elif Method == ComputationMethodsList[1]:  # Sodano
        Route = Geodesy.InverseSodano(OriginPoint=Origin,
                                      DestinationPoint=Dest)
    elif Method == ComputationMethodsList[2]:  #Sphere
        Route = Geodesy.InverseShperical(OriginPoint=Origin,
                                         DestinationPoint=Dest)
    else:
        messagebox.showerror(title="Invalid method",
                             message="Method " + Method + " not supported ATM")
        return
    TxtStartDist.insert(0, str("{:.2f}".format(Route.OrthoDistance / 1000)) + " km")
    TxtStartBear.insert(0, str("{:.2f}".format(math.degrees(Route.FwdAz))) + "°")
    TxtFinalAz.insert(0, str("{:.2f}".format(math.degrees(Route.BackAz))) + "°")


def Direct_Imp():
    Origin = Geodesy.Point()
    Dest = Geodesy.Point()
    Route = Geodesy.Route()
    try:
        Origin.Latitude = math.radians(float(TxtStartLat.get()))
        Origin.Longitude = math.radians(float(TxtStartLon.get()))
        Route.FwdAz = Geodesy.ParseAngle(TxtStartBear.get())
        Route.OrthoDistance = Geodesy.ParseDistance(TxtStartDist.get())
        if Route.OrthoDistance != Route.OrthoDistance or Route.FwdAz != Route.FwdAz:
            raise ValueError
    except ValueError:
        messagebox.showerror(
            title="Direct distance",
            message="Check origin point and route information"
        )
        return
    Method = ListBoxMethod.get()
    if Method == ComputationMethodsList[0]:  # Vincenty
        Dest = Geodesy.DirectVincenty(OriginPoint=Origin, Route=Route, tol=1e-24)
    elif Method == ComputationMethodsList[1]: #Sodano
        Dest = Geodesy.DirectSodano(OriginPoint=Origin, Route=Route)
    elif Method == ComputationMethodsList[2]: #Sphere
        Dest = Geodesy.DirectShperical(OriginPoint=Origin, Route=Route)
    else:
        messagebox.showerror(
            title="Invalid method", message="Method " + Method + " not supported ATM"
        )
        return
    TxtDestLat.delete(0, tkinter.END)
    TxtDestLat.insert(0, "{:.6f}".format((math.degrees(Dest.Latitude))) + " °")
    TxtDestLon.delete(0, tkinter.END)
    TxtDestLon.insert(0, "{:.6f}".format(math.degrees(Dest.Longitude))  + " °")


def SearchPoint():
    SearchKey: str = TxtSearch1.get()
    global SearchResult
    SearchResult.clear()
    SearchResult = MyDatabase.Search(key=SearchKey)
    if len(SearchResult) < 1:
        TxtStartLon.delete(0, tkinter.END)
        TxtStartLat.delete(0, tkinter.END)
        TxtStartLat.insert(0, "NOT FOUND")
        return
    elif len(SearchResult) == 1:
        TxtStartLat.delete(0, tkinter.END)
        TxtStartLat.insert(0, str(SearchResult[0].Lat))
        TxtStartLon.delete(0, tkinter.END)
        TxtStartLon.insert(0, str(SearchResult[0].Lon))
        return
    else:
        global MultipleSel
        MultipleSel = tkinter.Toplevel(master=home)
        MultipleSel.grab_set()
        global SearchResultBox
        SearchResultBox = tkinter.Listbox(
            master=MultipleSel, selectmode="single", font=DefaultFontTuple
        )
        for item in SearchResult:
            SearchResultBox.insert(tkinter.END, item.ICAO + "|" + item.Region)
        SearchResultBox.grid(row=0, column=0)
        global CmdSelect
        CmdSelect = tkinter.Button(
            master=MultipleSel,
            text="Select",
            font=DefaultFontTuple,
            command=ConfirmSelection,
        )
        CmdSelect.grid(row=1, column=0)
        return


def SearchPoint2():
    SearchKey: str = TxtSearch2.get()
    global SearchResult2
    SearchResult2.clear()
    SearchResult2 = MyDatabase.Search(key=SearchKey)
    if len(SearchResult2) < 1:
        TxtDestLon.delete(0, tkinter.END)
        TxtDestLat.delete(0, tkinter.END)
        TxtDestLat.insert(0, "NOT FOUND")
        return
    elif len(SearchResult2) == 1:
        TxtDestLat.delete(0, tkinter.END)
        TxtDestLat.insert(0, str(SearchResult2[0].Lat))
        TxtDestLon.delete(0, tkinter.END)
        TxtDestLon.insert(0, str(SearchResult2[0].Lon))
    else:
        global MultipleSel
        MultipleSel = tkinter.Toplevel(master=home)
        MultipleSel.grab_set()
        global SearchResultBox
        SearchResultBox = tkinter.Listbox(
            master=MultipleSel, selectmode="single", font=DefaultFontTuple
        )
        for item in SearchResult2:
            SearchResultBox.insert(tkinter.END, item.ICAO + "|" + item.Region)
        SearchResultBox.grid(row=0, column=0)
        global CmdSelect
        CmdSelect = tkinter.Button(
            master=MultipleSel,
            text="Select",
            font=DefaultFontTuple,
            command=ConfirmSelection2,
        )
        CmdSelect.grid(row=1, column=0)
        return
    return


home = tkinter.Tk(className="Geodesy calculator")
home.title("Geodesy calculator")
home.resizable(width=False, height=False)

InverseIcon = tkinter.PhotoImage(file="./icons/VinInverse.png")
InverseIcon = InverseIcon.subsample(x=3, y=3)
DirectIcon = tkinter.PhotoImage(file="./icons/VinDirect.png")
DirectIcon = DirectIcon.subsample(x=3, y=3)
SearchIcon = tkinter.PhotoImage(file="./icons/search.png")
SearchIcon = SearchIcon.subsample(x=3, y=3)
LoadDbIcon = tkinter.PhotoImage(file="./icons/database.png")
LoadDbIcon = LoadDbIcon.subsample(x=3, y=3)

VhfCount = tkinter.IntVar(value=0)
AptCount = tkinter.IntVar(value=0)
NdbCount = tkinter.IntVar(value=0)
Dbname = tkinter.StringVar(value="Database Name")

FrameDbInfo = tkinter.LabelFrame(
    master=home, text=Dbname.get(), font=DefaultFontTuple)
FrameDbInfo.grid(row=0, column=0, columnspan=2)
LblDbVhfCount = tkinter.Label(
    master=FrameDbInfo,
    text="VHF: " + str(VhfCount.get()).ljust(5, " "),
    font=DefaultFontTuple,
    foreground="green",
)
LblDbVhfCount.grid(row=0, column=0)
LblDbNdbCount = tkinter.Label(
    master=FrameDbInfo,
    text="NDB: " + str(NdbCount.get()).ljust(5, " "),
    font=DefaultFontTuple,
    foreground="red",
)
LblDbNdbCount.grid(row=0, column=1)
LblDbAptCount = tkinter.Label(
    master=FrameDbInfo,
    text="APT: " + str(AptCount.get()).ljust(5, " "),
    font=DefaultFontTuple,
    foreground="blue",
)
LblDbAptCount.grid(row=0, column=2)

CmdLoadDb = tkinter.Button(
    master=home,
    text="Load",
    font=DefaultFontTuple,
    command=LoadDb,
    image=LoadDbIcon,
    compound=tkinter.LEFT,
)
CmdLoadDb.grid(row=0, column=3)

FrameStart = tkinter.LabelFrame(
    master=home, text="Starting Point", font=DefaultFontTuple
)
FrameStart.grid(row=1, column=0, rowspan=2, padx=padValue, pady=padValue)
LblStartLat = tkinter.Label(
    master=FrameStart, text="Latitude", font=DefaultFontTuple)
LblStartLat.grid(row=0, column=0, columnspan=2)
TxtStartLat = tkinter.Entry(master=FrameStart, width=20, font=NumericFontTuple)
TxtStartLat.grid(row=1, column=0, columnspan=2)
LblStartLon = tkinter.Label(
    master=FrameStart, text="Longitude", font=DefaultFontTuple)
LblStartLon.grid(row=2, column=0, columnspan=2)
TxtStartLon = tkinter.Entry(master=FrameStart, width=20, font=NumericFontTuple)
TxtStartLon.grid(row=3, column=0, columnspan=2)

FrameRoute = tkinter.LabelFrame(
    master=home, text="Orthodromic route", font=DefaultFontTuple
)
FrameRoute.grid(row=1, column=1, columnspan=2, padx=padValue, pady=padValue)
LblStartBear = tkinter.Label(
    master=FrameRoute, text="Fwd. Az.", font=DefaultFontTuple)
LblStartBear.grid(row=0, column=0, columnspan=1)
TxtStartBear = tkinter.Entry(
    master=FrameRoute, width=CentralWidth, font=NumericFontTuple
)
TxtStartBear.grid(row=1, column=0, columnspan=1)
LblFinalAz = tkinter.Label(
    master=FrameRoute, text="Fin. Az.", font=DefaultFontTuple)
LblFinalAz.grid(row=0, column=1, columnspan=1)
TxtFinalAz = tkinter.Entry(
    master=FrameRoute, width=CentralWidth, font=NumericFontTuple)
TxtFinalAz.grid(row=1, column=1, columnspan=1)
LblStartDist = tkinter.Label(
    master=FrameRoute, text="Distance", font=DefaultFontTuple)
LblStartDist.grid(row=2, column=0, columnspan=1)
TxtStartDist = tkinter.Entry(
    master=FrameRoute, width=CentralWidth, font=NumericFontTuple
)
TxtStartDist.grid(row=3, column=0, columnspan=1)
LblMethod = tkinter.Label(
    master=FrameRoute, text="Method", font=DefaultFontTuple)
LblMethod.grid(row=2, column=1)
ListBoxMethod = ttk.Combobox(
    master=FrameRoute,
    font=DefaultFontTuple,
    state="readonly",
    width=CentralWidth,
    justify="center",
)
ListBoxMethod["values"] = ComputationMethodsList
ListBoxMethod.current(0)
ListBoxMethod.grid(row=3, column=1)

FrameDest = tkinter.LabelFrame(
    master=home, text="Destination Point", font=DefaultFontTuple
)
FrameDest.grid(row=1, column=3, rowspan=2, columnspan=2,
               padx=padValue, pady=padValue)
LblDestLat = tkinter.Label(
    master=FrameDest, text="Latitude", font=DefaultFontTuple)
LblDestLat.grid(row=0, column=0, columnspan=2)
TxtDestLat = tkinter.Entry(master=FrameDest, width=20, font=NumericFontTuple)
TxtDestLat.grid(row=1, column=0, columnspan=2)
LblDestLon = tkinter.Label(
    master=FrameDest, text="Longitude", font=DefaultFontTuple)
LblDestLon.grid(row=2, column=0, columnspan=2)
TxtDestLon = tkinter.Entry(master=FrameDest, width=20, font=NumericFontTuple)
TxtDestLon.grid(row=3, column=0, columnspan=2)

TxtSearch1 = tkinter.Entry(master=FrameStart, width=6, font=DefaultFontTuple)
TxtSearch1.grid(row=4, column=0)
CmdSearch1 = tkinter.Button(
    master=FrameStart,
    text="Search",
    font=DefaultFontTuple,
    command=SearchPoint,
    image=SearchIcon,
    compound=tkinter.LEFT,
    state="disabled",
)
CmdSearch1.grid(row=4, column=1)
CmdDirect = tkinter.Button(
    master=FrameRoute,
    text="Dir",
    font=DefaultFontTuple,
    image=DirectIcon,
    compound=tkinter.LEFT,
    command=Direct_Imp,
    state="active",
)
CmdDirect.grid(row=4, column=0)
CmdInverse = tkinter.Button(
    master=FrameRoute,
    text="Inv",
    font=DefaultFontTuple,
    command=Inverse_Imp,
    image=InverseIcon,
    compound=tkinter.LEFT,
)
CmdInverse.grid(row=4, column=1)
CmdSearch2 = tkinter.Button(
    master=FrameDest,
    text="Search",
    font=DefaultFontTuple,
    command=SearchPoint2,
    image=SearchIcon,
    compound=tkinter.LEFT,
    state="disabled",
)
CmdSearch2.grid(row=4, column=1)
TxtSearch2 = tkinter.Entry(master=FrameDest, width=6, font=DefaultFontTuple)
TxtSearch2.grid(row=4, column=0)


home.mainloop()
