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

    def __init__(self, Name: str = "", Lat: float = 0.0, Lon: float = 0.0) -> None:
        self.Latitude = np.float128(Lat)
        self.Longitude = np.float128(Lon)
        self.Name = Name


class Route:

    def __init__(self) -> None:
        self.FwdAz = np.float128(0.0)
        self.BackAz = np.float128(0.0)
        self.OrthoDistance = np.float128(0.0)


SemiMajorAxis = np.float128(6378388.0)#6378137.0)
f = np.float128(.003367003367)#1 / 298.257223563)
SemiMinorAxis = np.float128(6356911.946)#(1 - f) * SemiMajorAxis
scnd_Ecc_sqrd = (np.square(SemiMajorAxis) / np.square(SemiMinorAxis)) - 1


def InverseSodano(OriginPoint: Point, DestinationPoint: Point) -> Route:
    output = Route()
    DeltaLon = DestinationPoint.Longitude - OriginPoint.Longitude
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
    S_bo += a * ((f + f2) * sin_phi - f2 / 2 * phi2 * np.power(np.sin(phi), -1))
    S_bo += m * (-(f+f2)/2*phi - (f+f2)/2*sin_phi*cos_phi + f2/2*phi2*np.power(np.tan(phi), -1))
    S_bo += np.square(a) * -f2/2*sin_phi*cos_phi
    S_bo += np.square(m) * (f2/16*phi + f2/16*sin_phi*cos_phi - f2/2*phi2*np.power(np.tan(phi), -1) - f2/8*sin_phi*np.power(cos_phi, 3))
    S_bo += a*m* (f2 / 2 * phi2 * np.power(sin_phi, -1) + f2 / 2 * sin_phi * np.square(cos_phi))
    k1 = phi * (f + f2)
    k1 += a * (-f2 / 2 * sin_phi - f2 * phi2 * np.power(sin_phi, -1))
    k1 += m * (-1.25*f2*phi + f2/4*sin_phi*cos_phi + f2*phi2*np.power(np.tan(phi), -1))
    lam = k1 * c + DeltaLon
    cot_a12 = (sin_beta2 * cos_beta1 - np.cos(lam) * sin_beta1 * cos_beta2) / (np.sin(lam) * cos_beta2)
    cot_a21 = (sin_beta2 * cos_beta1 * np.cos(lam) - sin_beta1 * cos_beta2) / (np.sin(lam) * cos_beta1)
    alpha12 = np.arctan(np.power(cot_a12, -1))
    alpha21 = np.arctan(np.power(cot_a21, -1))
    output.OrthoDistance = S_bo * SemiMinorAxis
    output.FwdAz = alpha12
    output.BackAz = alpha21
    return output


def DirectSodano(OriginPoint: Point, Route: Route) -> Point:
    output = Point()
    sin_alpha12 = np.sin(Route.FwdAz)
    print("sin alpha 1-2 = " + str(sin_alpha12))
    cos_alpha12 = np.cos(Route.FwdAz)
    e2 = (np.square(SemiMajorAxis) / np.square(SemiMinorAxis)) - 1
    print("e2 = " + str(e2))
    e4 = np.square(e2)
    if np.absolute(OriginPoint.Latitude) <= (np.pi / 4):
        beta1 = np.arctan(math.tan(OriginPoint.Latitude) * (1 - f))
    else:
        cot_B = np.power(np.tan(OriginPoint.Latitude), -1)
        cot_beta1 = cot_B / (1 - f)
        beta1 = np.arctan(np.power(cot_beta1, -1))
    sin_beta1 = np.sin(beta1)
    print("sin beta1 = " + str(sin_beta1))
    cos_beta1 = np.cos(beta1)
    print("cos beta1 = " + str(cos_beta1))
    cos_beta0 = cos_beta1 * sin_alpha12
    g = cos_beta1 * cos_alpha12
    print("g = " + str(g))
    m1 = (1 + e2/2*np.square(sin_beta1)) * (1-np.square(cos_beta0))
    print("m1 = " + str(m1))
    phi_S = Route.OrthoDistance / SemiMinorAxis
    print("phiS = " + str(phi_S))
    sin_phiS = np.sin(phi_S)
    cos_phiS = np.cos(phi_S)
    a1 = (1 + e2/2*np.square(sin_beta1)) * (np.square(sin_beta1)*np.cos(phi_S) + g*sin_beta1*np.sin(phi_S))
    phi_0 = phi_S
    phi_0 += a1 * (-e2/2 * np.sin(phi_S))
    phi_0 += m1 * (e2/4 * (-phi_S+np.sin(phi_S)*np.cos(phi_S)))
    phi_0 += np.square(a1) * (0.625*e4*sin_phiS*cos_phiS)
    phi_0 += np.square(m1) * (11/64*e4*phi_S - 13/64*e4*sin_phiS*cos_phiS - e4/8*phi_S*np.square(cos_phiS) + 5/32*e4*sin_phiS*np.power(cos_phiS, 3))
    phi_0 += a1*m1 * (3/8*e4*sin_phiS + e4/4*phi_S*cos_phiS - 5/8*e4*sin_phiS*np.square(cos_phiS))
    sin_phi0 = np.sin(phi_0)
    cos_phi0 = np.cos(phi_0)
    cot_alpha21 = (g*cos_phi0 - sin_beta1*sin_phi0) / cos_beta0
    if np.absolute(sin_alpha12) < np.rad2deg(1e-6):
        alpha21 = 0
    else:
        if np.absolute(cot_alpha21) > 1:
            alpha21 = np.arctan(cot_alpha21)
        else:
            tan_alpha21 = 1 / cot_alpha21
            alpha21 = np.arctan(tan_alpha21)
    cot_lambda = (cos_beta1*cos_phi0 - sin_beta1*sin_phi0*cos_alpha12) / (sin_phi0*sin_alpha12)
    if np.absolute(sin_alpha12) < np.rad2deg(1e-6):
        lambda_ = 0
    else:
        if np.absolute(cot_lambda) > 1:
            lambda_ = np.arctan(cot_lambda)
        else:
            tan_lambda = 1 / cot_lambda
            lambda_ = np.arctan(tan_lambda)
    k = -f*phi_S + a1*1.5*f*f*sin_phiS + m1*(.75*f*f*phi_S - .75*f*f*sin_phiS*cos_phiS)
    L = k * cos_beta0 + lambda_
    print("L = " + str(L))
    output.Longitude = OriginPoint.Longitude + L
    if np.absolute(output.Longitude) > np.deg2rad(180.0):
        output.Longitude += -(2*np.pi) * np.sign(output.Longitude)
    sin_beta2 = sin_beta1*cos_phi0 + g*sin_phi0
    cos_beta2 = np.sqrt(np.square(cos_beta0) + np.square(g*cos_phi0-sin_beta1+sin_phi0))
    tan_beta2 = sin_beta2 / cos_beta2
    cot_beta2 = cos_beta2 / sin_beta2
    if np.absolute(tan_beta2) < np.absolute(cot_beta2):
        tan_B = tan_beta2 / (1-f)
        output.Latitude = np.arctan(tan_B)
    else:
        cot_B = cot_beta2 * (1-f)
        output.Latitude = np.arctan(1/cot_B)
    return output


def InverseVincenty(OriginPoint: Point, DestinationPoint: Point, tol: float = 1e-12) -> Route:
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
    O = Point(Name="origine", Lat=np.deg2rad(20.0), Lon = 0)
    D = Point(Name="destinazione", Lat=np.deg2rad(45.0), Lon =np.deg2rad(106.0))
    Rotta = InverseSodano(O, D)
    print("Sodano Inversa (meters) = " + str(Rotta.OrthoDistance)+"\n\n")
    Rotta2 = Route()
    Rotta2.FwdAz = np.deg2rad(42.94168)
    Rotta2.OrthoDistance = 9649412.505
    D = DirectSodano(OriginPoint=O,Route=Rotta)
    print("D = " + str(np.rad2deg(D.Latitude)) + ", " + str(np.rad2deg(D.Longitude)))
    exit()
