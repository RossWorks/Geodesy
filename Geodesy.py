import math
import numpy as np

DistanceUnits2Meters: dict = {
    "m": 1.0,
    "km": 1000.0,
    "nm": 1852.0,
    "mi": 1609.0,
}

Angle2RadConverter: dict = {"°": math.radians(1), "rad": 1.0}


def ParseDistance(UserInput: str) -> float:
    SpaceIndex = UserInput.find(" ")
    Number = float(UserInput[0:SpaceIndex])
    Unit = DistanceUnits2Meters[UserInput[SpaceIndex + 1 :].lower()]
    return Number * Unit


def ParseAngle(UserInput: str) -> float:
    SpaceIndex = UserInput.find(" ")
    Number = float(UserInput[0:SpaceIndex])
    Unit = Angle2RadConverter[UserInput[SpaceIndex + 1 :].lower()]
    return Number * Unit


class Point:
    Latitude: float
    Longitude: float

    def __init__(self, Name: str = "", Lat: float = 0.0, Lon: float = 0.0) -> None:
        self.Latitude = Lat
        self.Longitude = Lon
        self.Name = Name


class Route:
    FwdAz: float
    OrthoDistance: float
    BackAz: float

    def __init__(self) -> None:
        self.FwdAz = 0.0
        self.BackAz = 0.0
        self.OrthoDistance = 0.0


SemiMajorAxis: float = 6378137.0
f: float = 1 / 298.257223563
SemiMinorAxis: float = (1 - f) * SemiMajorAxis
scnd_Ecc_sqrd = (math.pow(SemiMajorAxis, 2) / math.pow(SemiMinorAxis, 2)) - 1


def InverseSodano(OriginPoint: Point, DestinationPoint: Point) -> Route:
    output = Route()
    DeltaLon = DestinationPoint.Longitude - OriginPoint.Longitude
    print("L = " + str(DeltaLon))
    if DeltaLon > np.pi:
        DeltaLon += -2 * np.pi * np.sign(DeltaLon)
    if np.absolute(OriginPoint.Latitude) <= (np.pi / 4):
        beta1 = np.arctan(math.tan(OriginPoint.Latitude) * (1 - f))
    else:
        cot_B = np.power(np.tan(OriginPoint.Latitude), -1)
        cot_beta1 = cot_B / (1 - f)
        beta1 = np.arctan(np.power(cot_beta1, -1))
    if np.absolute(DestinationPoint.Latitude) <= (np.pi / 4):
        beta2 = np.arctan(math.tan(DestinationPoint.Latitude) * (1 - f))
    else:
        cot_B = np.power(np.tan(DestinationPoint.Latitude), -1)
        cot_beta2 = cot_B / (1 - f)
        beta2 = np.arctan(np.power(cot_beta2, -1))
    sin_beta1 = np.sin(beta1)
    cos_beta1 = np.cos(beta1)
    sin_beta2 = np.sin(beta2)
    cos_beta2 = np.cos(beta2)
    a = sin_beta1 * sin_beta2
    b = cos_beta1 * cos_beta2
    cos_phi = a + b * np.cos(DeltaLon)
    sin_phi = sin_beta2 * cos_beta1 - sin_beta1 * cos_beta2 * np.cos(DeltaLon)
    sin_phi = np.power(sin_phi, 2)
    sin_phi = sin_phi + np.power(np.sin(DeltaLon) * cos_beta2, 2)
    sin_phi = np.sqrt(sin_phi)
    if np.absolute(sin_phi) < np.absolute(cos_phi):
        phi = np.arcsin(sin_phi)
    else:
        phi = np.arccos(cos_phi)
    c = b * np.sin(DeltaLon) / np.sin(phi)
    m = 1 - np.square(c)
    f2 = np.square(f)
    phi2 = np.square(phi)
    S_bo = phi * (1 + f + f2)
    S_bo += a * ((f + f2) * sin_phi - f2 / 2 * phi2 + np.power(np.cos(phi), -1))
    S_bo += m * (-(f+f2)/2*phi - (f+f2)/2*sin_phi*cos_phi + f2/2*phi2*np.power(np.tan(phi), -1))
    S_bo += np.square(a) * -f2 / 2 * sin_phi * cos_phi
    S_bo += np.square(m) * (f2/16*phi + f2/16*sin_phi*cos_phi - f2/2*phi2*np.power(np.tan(phi), -1) - f2/8*sin_phi*np.power(cos_phi, 3))
    S_bo += a*m* (f2 / 2 * phi2 * np.power(np.cos(phi), -1) + f2 / 2 * sin_phi * np.square(cos_phi))
    k1 = phi * (f + f2)
    k1 += a * (-f2 / 2 * sin_phi - f2 * phi2 * np.power(sin_phi, -1))
    k1 += m * (-1.25*f2 + f2/4*sin_phi*cos_phi + f2*phi2*np.power(np.tan(phi), -1))
    lam = k1 * c + DeltaLon
    cot_a12 = (sin_beta2 * cos_beta1 - np.cos(lam) * sin_beta1 * cos_beta2) / (np.sin(lam) * cos_beta2)
    cot_a21 = (sin_beta2 * cos_beta1 * np.cos(lam) - sin_beta1 * cos_beta2) / (np.sin(lam) * cos_beta1)
    alpha12 = np.arctan(np.power(cot_a12, -1))
    alpha21 = np.arctan(np.power(cot_a21, -1))
    output.OrthoDistance = S_bo * SemiMinorAxis
    output.FwdAz = alpha12
    output.BackAz = alpha21
    return output


def DirectSodano(OriginPoint: Point, FwdAz: float, Distance: float) -> Point:
    output = Point()
    return output


def InverseVincenty(
    OriginPoint: Point, DestinationPoint: Point, tol: float = 1e-12
) -> Route:
    """
    implemented form Wikipedia page
    https://en.wikipedia.org/wiki/Vincenty's_formulae
    inputs and outputs are in MKS system
    """
    Counter: int = 0
    output = Route()
    return output


def DirectVincenty(OriginPoint: Point, Route: Route, tol: float = 1e-12) -> Point:
    output = Point()

    return output


if __name__ == "__main__":
    print("Self test for " + __file__)
    Origins = [Point(Name="TOP"), Point(Name="TOP")]
    exit()
