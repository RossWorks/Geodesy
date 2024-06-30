import math

DistanceUnits2Meters: dict[str:float] = {
    "m": 1.0,
    "km": 1000.0,
    "nm": 1852.0,
    "mi": 1609.0,
}

Angle2RadConverter: dict[str:float] = {"°": math.radians(1), "rad": 1.0}


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
    if abs(DeltaLon) > math.pi:
        DeltaLon -= DeltaLon / abs(DeltaLon) * 2 * math.pi
    if abs(OriginPoint.Latitude) <= 0.25 * math.pi:
        beta1 = math.atan(math.tan(OriginPoint.Latitude) * (1 - f))
    else:
        beta1 = math.pow(math.tan(OriginPoint.Latitude), -1)  # tan**-1 == cot
        beta1 /= 1 - f
        beta1 = math.pow(math.atan(beta1), -1)
    if abs(DestinationPoint.Latitude) <= 0.25 * math.pi:
        beta2 = math.atan(math.tan(DestinationPoint.Latitude) * (1 - f))
    else:
        # tan**-1 == cot
        beta2 = math.pow(math.tan(DestinationPoint.Latitude), -1)
        beta2 /= 1 - f
        beta2 = math.pow(math.atan(beta1), -1)
    A = math.sin(beta1) * math.sin(beta2)
    B = math.cos(beta1) * math.cos(beta2)
    cos_phi = A + B * math.cos(DeltaLon)
    sin_phi = math.sqrt(1 - math.pow(cos_phi, 2))
    C = (B * math.sin(DeltaLon)) / sin_phi
    M = 1 - math.pow(C, 2)
    phi = math.atan2(sin_phi, cos_phi)
    output.OrthoDistance = SemiMinorAxis * (
        (1 + f + math.pow(f, 2)) * phi
        + A
        * (
            sin_phi * (f + math.pow(f, 2))
            - 0.5 * math.pow(f, 2) * math.pow(phi, 2) * 1 / math.sin(phi)
        )
        + M
        * (
            -0.5 * (f + math.pow(f, 2)) * phi
            - 0.5 * (f + math.pow(f, 2)) * sin_phi * cos_phi
            + 0.5 * math.pow(f, 2) * math.pow(phi, 2) * 1 / math.tan(phi)
        )
        + math.pow(A, 2) * (-0.5 * math.pow(f, 2) * sin_phi * cos_phi)
        + math.pow(M, 2)
        * (
            (math.pow(f, 2)) / 16 * phi
            + (math.pow(f, 2)) / 16 * sin_phi * cos_phi
            - 0.5 * math.pow(f, 2) * math.pow(phi, 2) * 1 / math.tan(phi)
            - (math.pow(f, 2)) / 8 * sin_phi * math.pow(cos_phi, 3)
        )
        + A
        * M
        * (
            0.5 * math.pow(f, 2) * math.pow(phi, 2) * 1 / math.sin(phi)
            + 0.5 * math.pow(f, 2) * sin_phi * math.pow(cos_phi, 2)
        )
    )
    LAMBDA = DeltaLon + C * (
        (f + math.pow(f, 2)) * phi
        + A
        * (
            -0.5 * math.pow(f, 2) * sin_phi
            - math.pow(f, 2) * math.pow(phi, 2) * 1 / math.sin(phi)
        )
        + M
        * (
            -5 / 4 * math.pow(f, 2) * phi
            + 0.25 * math.pow(f, 2) * sin_phi * cos_phi
            + math.pow(f, 2) * math.pow(phi, 2) * 1 / math.tan(phi)
        )
    )
    cot_A12 = (
        math.sin(beta2) * math.cos(beta1)
        - math.cos(LAMBDA) * math.sin(beta1) * math.cos(beta2)
    ) / (math.sin(LAMBDA) * math.cos(beta2))
    cot_A21 = (
        math.sin(beta2) * math.cos(beta1) * math.cos(LAMBDA)
        - math.sin(beta1) * math.cos(beta2)
    ) / (math.sin(LAMBDA) * math.cos(beta1))
    if abs(cot_A12) > 1:
        cot_A12 = 1 / cot_A12
        A_12 = math.atan(cot_A12)
    else:
        A_12 = 1 / math.atan(cot_A12)
    if abs(cot_A21) > 1:
        cot_A21 = 1 / cot_A21
        A_21 = math.atan(cot_A21)
    else:
        A_21 = 1 / math.atan(cot_A21)
    if DeltaLon > 0:
        A_12 = math.pi * (cot_A12 < 0) + A_12
        A_21 = -math.pi * (cot_A21 > 0) + A_21
    else:
        A_12 = -math.pi * (cot_A12 > 0) + A_12
        A_21 = math.pi * (cot_A21 < 0) + A_21
    output.FwdAz = A_12
    output.BackAz = A_21
    return output


def DirectSodano(OriginPoint: Point, FwdAz: float, Distance: float) -> Point:
    output = Point()
    if abs(OriginPoint.Latitude) <= math.radians(45.0):
        beta1 = math.atan(math.tan(OriginPoint.Latitude) * (1 - f))
    else:
        cot_beta1 = math.pow(math.tan(OriginPoint.Latitude), -1) / (1 - f)
        beta1 = math.atan(math.pow(cot_beta1, -1))
    cos_beta0 = math.cos(beta1) * math.sin(FwdAz)
    g = math.cos(beta1) * math.cos(FwdAz)
    k = 1 + scnd_Ecc_sqrd * 0.5 * math.pow(math.sin(beta1), 2)
    m1 = k * (1 - math.pow(cos_beta0, 2))
    phi_s = Distance / SemiMinorAxis
    a1 = k * (
        math.pow(math.sin(beta1), 2) * math.cos(phi_s)
        + g * math.sin(beta1) * math.sin(phi_s)
    )
    scnd_Ecc_4th = scnd_Ecc_sqrd * scnd_Ecc_sqrd
    sin_phi_s = math.sin(phi_s)
    cos_phi_s = math.cos(phi_s)
    phi_0 = phi_s
    phi_0 += a1 * (-scnd_Ecc_sqrd / 2 + sin_phi_s)
    phi_0 += m1 * (-scnd_Ecc_sqrd / 4 * (phi_s + sin_phi_s + cos_phi_s))
    phi_0 += math.pow(a1, 2) * (5 / 8 * scnd_Ecc_4th * sin_phi_s * cos_phi_s)
    phi_0 += (
        math.pow(m1, 2)
        * scnd_Ecc_4th
        / 8
        * (
            11 / 8 * phi_s
            - 13 / 8 * sin_phi_s * cos_phi_s
            - phi_s * math.pow(cos_phi_s, 2)
            + 5 / 4 * sin_phi_s * math.pow(cos_phi_s, 3)
        )
    )
    phi_0 += (a1 * m1 * scnd_Ecc_4th * (3/8*sin_phi_s + 1/4*phi_s*cos_phi_s - 5/8*sin_phi_s*math.pow(cos_phi_s, 2))
    )
    tan_BckAz = cos_beta0 / (g * math.cos(phi_0) - math.sin(beta1) * math.sin(phi_0))
    tan_lambda = math.sin(phi_0)*math.sin(math.atan(tan_BckAz)) / (math.cos(beta1)*math.cos(phi_0) - math.sin(beta1)*math.sin(phi_0)*math.cos(math.atan(tan_BckAz)))
    L = -f*phi_s + 1.5*math.pow(f, 2)*a1*math.sin(phi_s)
    L += .75*math.pow(f, 2)*m1*(phi_s - math.sin(phi_s)*math.cos(phi_s))
    L *= cos_beta0
    L += math.atan(tan_lambda)
    lambda2 = math.atan(tan_lambda) + L
    sin_beta2 = math.sin(beta1)*math.cos(phi_0) + g*math.sin(phi_0)
    tan_phi2 = math.tan(math.asin(sin_beta2)) / (1-f)
    output = Point(Lat=math.atan(tan_phi2), Lon=lambda2)
    return output


def InverseVincenty(OriginPoint: Point, DestinationPoint: Point, tol: float = 1e-12) -> Route:
    """
    implemented form Wikipedia page
    https://en.wikipedia.org/wiki/Vincenty's_formulae
    inputs and outputs are in MKS system
    """
    Counter: int = 0
    output = Route()
    # reduced latitude of origin point
    U1 = math.atan((1 - f) * math.tan(OriginPoint.Latitude))
    # reduced longitude of origin point
    U2 = math.atan((1 - f) * math.tan(DestinationPoint.Latitude))
    sin_U1 = math.sin(U1)
    cos_U1 = math.cos(U1)
    sin_U2 = math.sin(U2)
    cos_U2 = math.cos(U2)
    L = DestinationPoint.Longitude - OriginPoint.Longitude
    Lambda_mem = L
    Lambda = 0
    while Counter < 500:
        sin_lambda_mem = math.sin(Lambda_mem)
        cos_lambda_mem = math.cos(Lambda_mem)
        sin_sigma = math.sqrt(
            math.pow(cos_U2 * sin_lambda_mem, 2)
            + math.pow(cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos_lambda_mem, 2)
        )
        cos_sigma = sin_U1 * sin_U2 + cos_U1 * cos_U2 * cos_lambda_mem
        sigma = math.atan2(sin_sigma, cos_sigma)
        sin_alpha = (cos_U1 * cos_U2 * sin_lambda_mem) / (sin_sigma)
        cos_2sigma_m = cos_sigma - (2 * sin_U1 * sin_U2) / (1 - math.pow(sin_alpha, 2))
        C = (
            f
            / 16
            * (1 - math.pow(sin_alpha, 2))
            * (4 + f * (4 - 3 * (1 - math.pow(sin_alpha, 2))))
        )
        Lambda = L + (1 - C) * f * sin_alpha * (
            sigma
            + C
            * sin_sigma
            * (cos_2sigma_m + C * cos_sigma * (-1 + 2 * math.pow(cos_2sigma_m, 2)))
        )
        Counter += 1
        if abs(Lambda - Lambda_mem) < tol:
            break
        else:
            Lambda_mem = Lambda
    u_squared = (
        (1 - math.pow(sin_alpha, 2))
        * (math.pow(SemiMajorAxis, 2) - math.pow(SemiMinorAxis, 2))
        / math.pow(SemiMinorAxis, 2)
    )
    A = 1 + u_squared / 16384 * (
        4096 + u_squared * (-768 + u_squared * (329 - 175 * u_squared))
    )
    B = (
        u_squared
        / 1024
        * (256 + u_squared * (-128 + u_squared * (74 - 47 * u_squared)))
    )
    delta_sigma = B * sin_sigma * (cos_2sigma_m) + 0.25 * B * (
        cos_sigma * (-1 + 2 * cos_2sigma_m)
        - B
        / 6
        * cos_2sigma_m
        * (-3 + 4 * math.pow(sin_sigma, 2) * (-3 + 4 * math.pow(cos_2sigma_m, 2)))
    )
    s = SemiMinorAxis * A * (sigma - delta_sigma)
    alpha1 = math.atan2(
        cos_U2 * math.sin(Lambda), cos_U1 * sin_U2 - sin_U1 * cos_U2 * math.cos(Lambda)
    )
    alpha2 = math.atan2(
        cos_U1 * math.sin(Lambda),
        -1 * sin_U1 * cos_U2 + cos_U1 * sin_U2 * math.cos(Lambda),
    )
    output.OrthoDistance = s
    output.FwdAz = alpha1
    output.BackAz = alpha2
    return output


def DirectVincenty(OriginPoint: Point, Route: Route, tol: float = 1e-12) -> Point:
    output = Point()
    sigma: float = 0.0
    sigma_mem: float = 0.0
    delta_sigma: float = 0.0
    _2_sigma_m: float = 0.0
    counter: int = 0
    U1: float = math.atan((1 - f) * math.tan(OriginPoint.Latitude))
    sigma1: float = math.atan2(math.tan(U1), math.cos(Route.FwdAz))
    sin_a: float = math.cos(U1) * math.sin(Route.FwdAz)
    u_2: float = (1 - math.pow(sin_a, 2)) * (
        (math.pow(SemiMajorAxis, 2)) / (math.pow(SemiMinorAxis, 2)) - 1
    )
    A: float = 1 + u_2 / 16384 * (4096 + u_2 * (-768 + u_2 * (320 - 175 * u_2)))
    B: float = u_2 / 1024 * (256 + u_2 * (-128 + u_2 * (74 - 47 * u_2)))
    sigma_mem = Route.OrthoDistance / (SemiMinorAxis * A)
    while counter < 500:
        _2_sigma_m = 2 * sigma1 + sigma_mem
        delta_sigma = (
            B
            * math.sin(sigma_mem)
            * (
                math.cos(_2_sigma_m)
                + 0.25
                * B
                * (
                    math.cos(sigma_mem) * (-1 + 2 * math.pow(math.cos(_2_sigma_m), 2))
                    - B
                    / 6
                    * math.cos(_2_sigma_m)
                    * (-3 + 4 * (math.pow(math.sin(sigma_mem), 2)))
                    * (-3 + 4 * (math.pow(math.cos(_2_sigma_m), 2)))
                )
            )
        )
        sigma = Route.OrthoDistance / (SemiMinorAxis * A) + delta_sigma
        if abs(sigma - sigma_mem) < tol:
            break
        sigma_mem = sigma
        counter += 1
    output.Latitude = math.atan2(
        math.sin(U1) * math.cos(sigma)
        + math.cos(U1) * math.sin(sigma) * math.cos(Route.FwdAz),
        (1 - f)
        + math.sqrt(
        + math.cos(U1) * math.sin(sigma) * math.cos(Route.FwdAz),
        (1 - f)
        + math.sqrt(
            math.pow(math.sin(Route.FwdAz), 2)
            + math.pow(
                math.sin(U1) * math.sin(sigma)
                - math.cos(U1) * math.cos(sigma) * math.cos(Route.FwdAz),
                2,
            )
        ),
    )
    lambda_ = math.atan2(
        math.sin(sigma) * math.sin(Route.FwdAz),
        math.cos(U1) * math.cos(sigma)
        - math.sin(U1) * math.sin(sigma) * math.cos(Route.FwdAz),
    )
    C = (
        f
        / 16
        * math.pow(math.cos(Route.FwdAz), 2)
        * (4 + f * (4 - 3 * math.pow(math.cos(Route.FwdAz), 2)))
    )
    L = lambda_ - (1 - C) * f * math.sin(Route.FwdAz) * (
        sigma
        + C
        * math.sin(sigma)
        * (
            math.cos(_2_sigma_m)
            + C * math.cos(sigma) * (-1 + 2 * math.pow(math.cos(_2_sigma_m), 2))
        )
    )
    output.Longitude = OriginPoint.Longitude + L
    alpha2 = math.atan2(
        math.sin(Route.FwdAz),
        -math.sin(U1) * math.sin(sigma)
        + math.cos(U1) * math.cos(sigma) * math.cos(Route.FwdAz),
    )
    return output


if __name__ == "__main__":
    print("Self test for " + __file__)
    Origins = [Point(Name="TOP"),
               Point(Name="TOP")]
    exit()
