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
  def __init__(self) -> None:
    self.Latitude = np.float128(0.0)
    self.Longitude = np.float128(0.0)

class Route:
  def __init__(self) -> None:
    self.FwdAz = np.float128(0.0)
    self.BackAz = np.float128(0.0)
    self.OrthoDistance = np.float128(0.0)

SemiMajorAxis : np.float128 = 6378137.0
flattening    : np.float128 = 1/298.257223563
SemiMinorAxis : np.float128 = (1-flattening)*SemiMajorAxis

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
  L = DestinationPoint.Longitude-OriginPoint.Longitude
  Lambda_mem = L
  Lambda = 0
  while Counter < 500:
    sin_lambda_mem = np.sin(Lambda_mem)
    cos_lambda_mem = np.cos(Lambda_mem)
    sin_sigma = np.sqrt(np.power(cos_U2*sin_lambda_mem,2) + np.power(cos_U1*sin_U2 - sin_U1*cos_U2*cos_lambda_mem, 2))
    cos_sigma = sin_U1*sin_U2 + cos_U1*cos_U2*cos_lambda_mem
    sigma = np.arcsin(sin_sigma)
    sin_alpha = (cos_U1*cos_U2*sin_lambda_mem)/(sin_sigma)
    try:
      cos_2sigma_m = cos_sigma - (2*sin_U1*sin_U2)/(1-np.power(sin_alpha,2))
    except ZeroDivisionError:
      cos_2sigma_m = 0.0
    C = flattening/16 * (1-np.power(sin_alpha,2)) * (4+flattening*(4-3*(1-np.power(sin_alpha,2))))
    Lambda = L + (1-C) * flattening * sin_alpha * (sigma + C*sin_sigma * (cos_2sigma_m + C*cos_sigma*(-1+2*np.power(cos_2sigma_m,2))))
    Counter += 1
    if abs(Lambda - Lambda_mem) < tol:
      break
    else:
      Lambda_mem = Lambda
  u_squared = (1-np.power(sin_alpha,2))*(np.power(SemiMajorAxis,2)/np.power(SemiMinorAxis,2) - 1.0)
  A = 1 + u_squared / 16384 * (4096 + u_squared * (-768 + u_squared * (320 - 175 * u_squared)))
  B = u_squared / 1024 * (256 + u_squared * (-128 + u_squared * (74 - 47 * u_squared)))
  delta_sigma = B * sin_sigma * (cos_2sigma_m) + 0.25 * B * (cos_sigma * (-1 + 2 * np.power(cos_2sigma_m, 2)) - B / 6 *cos_2sigma_m*(-3+4*np.power(sin_sigma,2)*(-3+4*np.power(cos_2sigma_m,2))))
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