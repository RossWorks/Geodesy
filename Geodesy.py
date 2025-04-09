if __name__ == "__main__":
  print("Standalone operation not available")
  exit()

import numpy as np

DistanceUnits2Meters : dict[str:np.float128] = {'m'  : np.float128(    1.0),
                                                'km' : np.float128( 1000.0),
                                                'nm' : np.float128( 1852.0),
                                                'mi' : np.float128( 1609.0)}

Angle2RadConverter : dict[str:np.float128] = {"°"   : np.float128(np.radians(1)),
                                              "rad" : 1.0}

def ParseDistance(UserInput : str) -> np.float128:
  SpaceIndex = UserInput.find(' ')
  Number = np.float128(UserInput[0:SpaceIndex])
  Unit = DistanceUnits2Meters[UserInput[SpaceIndex+1:].lower()]
  return Number * Unit

def ParseAngle(UserInput : str) -> np.float128:
  SpaceIndex = UserInput.find(' ')
  Number = np.float128(UserInput[0:SpaceIndex])
  Unit = Angle2RadConverter[UserInput[SpaceIndex+1:].lower()]
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

# as per WGS84
SemiMajorAxis = np.float128(6378388.0)#6378137.0)
flattening = np.float128(.003367003367)#1 / 298.257223563)
SemiMinorAxis = np.float128(6356911.946)#(1 - f) * SemiMajorAxis
scnd_Ecc_sqrd = (np.square(SemiMajorAxis) / np.square(SemiMinorAxis)) - 1


def InverseSodano(OriginPoint: Point, DestinationPoint: Point) -> Route:
    output = Route()
    DeltaLon = DestinationPoint.Longitude - OriginPoint.Longitude
    if DeltaLon > np.pi:
        DeltaLon += -2 * np.pi * np.sign(DeltaLon)
    if np.absolute(OriginPoint.Latitude) <= (np.pi / 4):
        beta1 = np.arctan(np.tan(OriginPoint.Latitude) * (1 - flattening))
    else:
        cot_B = np.power(np.tan(OriginPoint.Latitude), -1)
        cot_beta1 = cot_B / (1 - flattening)
        beta1 = np.arctan(np.power(cot_beta1, -1))
    if np.absolute(DestinationPoint.Latitude) <= (np.pi / 4):
        beta2 = np.arctan(np.tan(DestinationPoint.Latitude) * (1 - flattening))
    else:
        cot_B = np.power(np.tan(DestinationPoint.Latitude), -1)
        cot_beta2 = cot_B / (1 - flattening)
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
    f2 = np.square(flattening)
    phi2 = np.square(phi)
    S_bo = phi * (1 + flattening + f2)
    S_bo += a * ((flattening + f2) * sin_phi - f2 / 2 * phi2 * np.power(np.sin(phi), -1))
    S_bo += m * (-(flattening+f2)/2*phi - (flattening+f2)/2*sin_phi*cos_phi + f2/2*phi2*np.power(np.tan(phi), -1))
    S_bo += np.square(a) * -f2/2*sin_phi*cos_phi
    S_bo += np.square(m) * (f2/16*phi + f2/16*sin_phi*cos_phi - f2/2*phi2*np.power(np.tan(phi), -1) - f2/8*sin_phi*np.power(cos_phi, 3))
    S_bo += a*m* (f2 / 2 * phi2 * np.power(sin_phi, -1) + f2 / 2 * sin_phi * np.square(cos_phi))
    k1 = phi * (flattening + f2)
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
        beta1 = np.arctan(np.tan(OriginPoint.Latitude) * (1 - flattening))
    else:
        cot_B = np.power(np.tan(OriginPoint.Latitude), -1)
        cot_beta1 = cot_B / (1 - flattening)
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
    k = -flattening*phi_S + a1*1.5*flattening*flattening*sin_phiS + m1*(.75*flattening*flattening*phi_S - .75*flattening*flattening*sin_phiS*cos_phiS)
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
        tan_B = tan_beta2 / (1-flattening)
        output.Latitude = np.arctan(tan_B)
    else:
        cot_B = cot_beta2 * (1-flattening)
        output.Latitude = np.arctan(1/cot_B)
    return output


def InverseVincenty(OriginPoint : Point, DestinationPoint : Point, tol : np.float128 = 1e-12) -> Route:
  '''
  implemented form Wikipedia page
  https://en.wikipedia.org/wiki/Vincenty's_formulae
  inputs and outputs are in MKS system
  '''
  Counter : int = 0
  output = Route()
  U1 = np.arctan((1-flattening)*np.tan(OriginPoint.Latitude)) #reduced latitude of origin point
  U2 = np.arctan((1-flattening)*np.tan(DestinationPoint.Latitude)) #reduced longitude of origin point
  sin_U1 = np.sin(U1)
  cos_U1 = np.cos(U1)
  sin_U2 = np.sin(U2)
  cos_U2 = np.cos(U2)
  L = DestinationPoint.Longitude - OriginPoint.Longitude
  Lambda_mem = L
  Lambda = 0
  for Counter in range(500):
    sin_lambda_mem = np.sin(Lambda_mem)
    cos_lambda_mem = np.cos(Lambda_mem)
    sin_sigma = np.sqrt(np.square(cos_U2*sin_lambda_mem) + np.square(cos_U1*sin_U2 - sin_U1*cos_U2*cos_lambda_mem))
    cos_sigma = sin_U1*sin_U2 + cos_U1*cos_U2*cos_lambda_mem
    sigma = np.arctan2(sin_sigma, cos_sigma)
    sin_alpha = (cos_U1*cos_U2*sin_lambda_mem)/(sin_sigma)
    cos2_alpha = 1 - np.square(sin_alpha)
    try:
      cos_2sigma_m = cos_sigma - (2*sin_U1*sin_U2)/cos2_alpha
    except ZeroDivisionError:
      cos_2sigma_m = 0.0
    C = flattening/16 * cos2_alpha * (4+flattening*(4-3*cos2_alpha))
    Lambda = L + (1-C) * flattening * sin_alpha * (sigma + C*sin_sigma * (cos_2sigma_m + C*cos_sigma*(-1+2*np.power(cos_2sigma_m,2))))
    if abs(Lambda - Lambda_mem) < tol:
      break
    else:
      Lambda_mem = Lambda
  u_squared = (1-np.power(sin_alpha,2))*(np.power(SemiMajorAxis,2)/np.power(SemiMinorAxis,2) - 1.0)
  A = 1 + u_squared / 16384 * (4096 + u_squared * (-768 + u_squared * (320 - 175 * u_squared)))
  B = u_squared / 1024 * (256 + u_squared * (-128 + u_squared * (74 - 47 * u_squared)))
  delta_sigma = B * sin_sigma * (cos_2sigma_m + 0.25 * B * (cos_sigma * (-1 + 2 * np.power(cos_2sigma_m, 2)) - B / 6 *cos_2sigma_m*(-3+4*np.power(sin_sigma,2)*(-3+4*np.power(cos_2sigma_m,2)))))
  s = SemiMinorAxis * A * (sigma-delta_sigma)
  alpha1 = np.arctan2(cos_U2 * np.sin(Lambda), cos_U1 * sin_U2 - sin_U1 * cos_U2 * np.cos(Lambda))
  alpha2 = np.arctan2(cos_U1 * np.sin(Lambda), -1*sin_U1*cos_U2+cos_U1*sin_U2*np.cos(Lambda))
  output.OrthoDistance = s
  output.FwdAz = alpha1
  output.BackAz = alpha2
  return output

def DirectVincenty(OriginPoint : Point, Route : Route, tol : np.float128 = 1e-12) -> Point:
  output = Point()
  sigma : np.float128 = 0.0
  sigma_mem : np.float128 = 0.0
  delta_sigma : np.float128 = 0.0
  _2_sigma_m : np.float128 = 0.0
  counter: int = 0
  U1     : np.float128 = np.arctan((1-flattening)*np.tan(OriginPoint.Latitude))
  sigma1 : np.float128 = np.arctan2(np.tan(U1), np.cos(Route.FwdAz))
  sin_a  : np.float128 = np.cos(U1)*np.sin(Route.FwdAz)
  u_2    : np.float128 = (1-np.power(sin_a,2))*((np.power(SemiMajorAxis,2))/(np.power(SemiMinorAxis,2))-1)
  A      : np.float128 = 1 + u_2/16384*(4096+u_2*(-768+u_2*(320-175*u_2)))
  B      : np.float128 = u_2/1024 * (256+u_2*(-128 + u_2*(74 - 47*u_2)))
  sigma_mem = Route.OrthoDistance/(SemiMinorAxis*A)
  while counter < 500:
    _2_sigma_m = 2*sigma1+sigma_mem
    delta_sigma = B * np.sin(sigma_mem)*(np.cos(_2_sigma_m) + .25*B * (np.cos(sigma_mem)*(-1+2*np.power(np.cos(_2_sigma_m),2)) - B/6*np.cos(_2_sigma_m) * (-3+4*(np.power(np.sin(sigma_mem),2))) * (-3+4*(np.power(np.cos(_2_sigma_m),2)))))
    sigma = Route.OrthoDistance/(SemiMinorAxis*A) + delta_sigma
    if (abs(sigma-sigma_mem) < tol):
      break
    sigma_mem = sigma
    counter += 1
  output.Latitude = np.arctan2(np.sin(U1)*np.cos(sigma) + np.cos(U1)*np.sin(sigma)*np.cos(Route.FwdAz),
                               (1-flattening)+np.sqrt(np.power(np.sin(Route.FwdAz),2) + np.power(np.sin(U1)*np.sin(sigma) - np.cos(U1)*np.cos(sigma)*np.cos(Route.FwdAz),2)))
  lambda_ = np.arctan2(np.sin(sigma)*np.sin(Route.FwdAz),
                       np.cos(U1)*np.cos(sigma) - np.sin(U1)*np.sin(sigma)*np.cos(Route.FwdAz))
  C = flattening/16 * np.power(np.cos(Route.FwdAz),2) * (4+flattening*(4-3*np.power(np.cos(Route.FwdAz),2)))
  L = lambda_ - (1-C)*flattening*np.sin(Route.FwdAz) * (sigma + C*np.sin(sigma)*(np.cos(_2_sigma_m) + C*np.cos(sigma)*(-1+2*np.power(np.cos(_2_sigma_m),2))))
  output.Longitude = OriginPoint.Longitude + L
  alpha2 = np.arctan2(np.sin(Route.FwdAz),
                      -np.sin(U1)*np.sin(sigma) + np.cos(U1)*np.cos(sigma)*np.cos(Route.FwdAz))
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