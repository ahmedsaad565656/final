import streamlit as st
import math
def calculate_specific_gravity():
    st.title("Oil Specific Gravity “γo” ")


    # Ask the user for the oil density
    rho_o = st.number_input("Enter the oil density (in lb/ft^3)", min_value=0.0)

    # Calculate the oil specific gravity and display the result
    if st.button("Calculate γo"):
        gamma_o = rho_o / 62.4
        st.success(f"Oil specific gravity: {gamma_o:.4f}")
#.API Gravity.
def calculate_API_gravity():
    st.title("“API” Gravity  ")


    # Ask the user for the specific gravity
    γo = st.number_input("Enter the oil specific gravity", min_value=0.0)

    # Calculate the oil specific gravity and display the result
    if st.button("Calculate API"):
        API = (141.5/γo)-131.5
        st.success(f"API gravity: {API:.4f}°API")
#..
def calculate_Gas_Solubility_McCain():
    st.title("Gas Solubility “Rs”using “McCain” Correlation")
    Bo = st.number_input("Enter the oil formation volume factor (bbl/STB)", min_value=0.0)
    ρo = st.number_input("Enter the oil density (lb/ft3)", min_value=0.0)
    Yo=st.number_input("Enter the oil specific gravity ", min_value=0.0)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.0)

    if st.button("Calculate Rs McCain"):

        Rs = (Bo*ρo-62.4*Yo)/(.0136*Yg)

        st.success(f"Gas solubility: {Rs:.4f}")


def calculate_Gas_Solubility_standing():
    st.title("Gas Solubility “Rs”using “Standing” Correlation")
    p = st.number_input("Enter the pressure (psia)", min_value=0.0)
    temperature = st.number_input("Enter the temperature (R)", min_value=0.0)
    API=st.number_input("Enter the API ", min_value=0.0)
    Yg = st.number_input("Enter the specific gravity ", min_value=0.0)

    if st.button("Calculate Rs standing"):
        x = 0.0125 * API - 0.00091 * (temperature - 460)
        Rs = (Yg * (p / 18.2 + 1.4) * 10 ** x) ** 1.2048

        st.success(f"Gas solubility: {Rs:.4f}")

def calculate_gas_solubility_Vasquez():
    st.title("Gas Solubility “Rs” using “Vasquez-Beggs” Correlation ")

    API = st.number_input("Enter the oil API gravity", min_value=0.0)
    Tsep = st.number_input("Enter the separator temperature (in °R)", min_value=0.0)
    Psep = st.number_input("Enter the separator pressure (in psi)", min_value=0.0)
    gamma_g = st.number_input("Enter the specific gravity of the gas", min_value=0.0)
    pressure = st.number_input("Enter the pressure (in psi)", min_value=0.0)
    temperature = st.number_input("Enter the temperature (in °R)", min_value=0.0)

    if st.button("Calculate Rs vasquez"):
                C1 = 0.0362 if API <= 30 else 0.0178
                C2 = 1.0937 if API <= 30 else 1.1870
                C3 = 25.7240 if API <= 30 else 23.931

                gamma_gs = gamma_g * (1 + (5.912 * 10 ** -5) * API * (Tsep - 460) * math.log10(Psep / 114.7))
                Rs = C1 * gamma_gs * pressure ** C2 * math.exp(C3 * (API / temperature))

                st.success(f"Gas solubility: {Rs:.4f} scf/STB")
def calculate_Gas_Solubility_Marhoun():
    st.title("Gas Solubility “Rs” using “Marhoun” Correlation")

    p = st.number_input("Enter the pressure (psia)", min_value=0.2)
    T = st.number_input("Enter the temperature (R)", min_value=0.2)
    γg = st.number_input("Enter the specific gravity of gas", min_value=0.2)
    γo = st.number_input("Enter the specific gravity of oil", min_value=0.2)

    if st.button("Calculate Rs Marhoun"):
        a = 185.843208
        b = 1.877840
        c = -3.1437
        d = -1.32657
        e = 1.398441
        Rs = (a * γg ** b * γo ** c * T ** d * p) ** e
        st.success(f"Gas solubility: {Rs:.4f} scf/STB")
#...
def calculate_Gas_Solubility_Petrosky():
    st.title("Gas Solubility “Rs” using “Petrosky” Correlation")
    p = st.number_input("Enter the pressure (psia)", min_value=0.3)
    temperature = st.number_input("Enter the temperature (R)", min_value=0.3)
    API=st.number_input("Enter the API ", min_value=0.3)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.3)

    if st.button("Calculate Rs Petrosky"):
        x = 7.916*(10**-4)*API**1.5410 - 4.561*(10**-5)*(temperature - 460)**1.3911
        Rs =(((p/112.727)+12.340)*(Yg**.8439)*10**x)**1.73184

        st.success(f"Gas solubility Petrosky : {Rs:.4f}")
#...
def calculate_Gas_Solubility_Glaso():
    st.title("Gas Solubility “Rs” using “Glaso” Correlation")
    p = st.number_input("Enter the pressure (psia)", min_value=0.3)
    T = st.number_input("Enter the temperature (R)", min_value=0.3)
    API=st.number_input("Enter the API ", min_value=0.3)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.3)

    if st.button("Calculate Rs Glaso"):
        x = 2.8869-(14.1811-3.3093*math.log(p,10))**.5
        Pb_star=10**x
        Rs =Yg*(((API**.989)/(T-460)**.172)*Pb_star)**1.2255

        st.success(f"Gas solubility Glaso : {Rs:.4f}")
#...
def calculate_Bubble_point_pressure_Standing():
    st.title("Bubble point pressure “Pb” using “Standing” Correlation")
    Rs = st.number_input("Enter the Gas solubility ", min_value=0.4)
    T = st.number_input("Enter the temperature (R)", min_value=0.4)
    API=st.number_input("Enter the API ", min_value=0.4)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.4)

    if st.button("Calculate Pb Standing"):
        a=.00091*(T-460)-.0125*API
        Pb=18.2*((Rs/Yg)**.83*10**a-1.4)
        st.success(f"Bubble point pressure standing : {Pb:.4f}")
#...
def calculate_Bubble_point_pressure_Vasquez():
    st.title("Bubble point pressure “Pb” using “Vasquez” Correlation")
    Rs = st.number_input("Enter the Gas solubility ", min_value=0.4)
    T = st.number_input("Enter the temperature (R)", min_value=0.4)
    API=st.number_input("Enter the API ", min_value=0.4)
    Ygs = st.number_input("Enter the gas specific gravity at separator pressure ", min_value=0.4)

    if st.button("Calculate Pb Vasquez"):
        C1 = 27.624 if API <= 30 else 56.18
        C2 = .914328 if API <= 30 else .84246
        C3 = 11.172 if API <= 30 else 10.393
        a=-C3*API/T
        Pb=((C1*Rs/Ygs)*10**a)**C2
        st.success(f"Bubble point pressure Vasquez : {Pb:.4f}")
#...
def calculate_Bubble_point_pressure_Glaso():
    st.title("Bubble point pressure “Pb” using “Glaso” Correlation")
    Rs = st.number_input("Enter the Gas solubility ", min_value=0.4)
    t = st.number_input("Enter the temperature (F)", min_value=0.4)
    API=st.number_input("Enter the API ", min_value=0.4)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.4)

    if st.button("Calculate Pb Glaso"):
        a=.816
        b=.172
        c=-.989
        Pb_star=((Rs/Yg)**a)*(t**b)*(API**c)
        log_Pb=1.7669+1.7447*math.log10(Pb_star)-.30218*(Pb_star)**2
        Pb=10**log_Pb

        st.success(f"Bubble point pressure Glaso : {Pb:.4f}")
#..
def calculate_Bubble_point_pressure_Marhoun():
    st.title("Bubble point pressure “Pb” using “Marhoun” Correlation")
    Rs = st.number_input("Enter the Gas solubility ", min_value=0.4)
    T = st.number_input("Enter the temperature (R)", min_value=0.4)
    Yo=st.number_input("Enter the oil specific gravity ", min_value=0.4)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.4)

    if st.button("Calculate Pb Marhoun"):
        a=5.38088*10**-3
        b=.715082
        c=-1.87784
        d=3.1437
        e=1.32657
        Pb=a*(Rs**b)*(Yg**c)*(Yo**d)*(T**e)

        st.success(f"Bubble point pressure Marhoun : {Pb:.4f}")
#..
def calculate_Bubble_point_pressure_Petrosky():
    st.title("Bubble point pressure “Pb” using “Petrosky” Correlation")
    Rs = st.number_input("Enter the Gas solubility", min_value=0.3)
    T = st.number_input("Enter the temperature (R)", min_value=0.3)
    API=st.number_input("Enter the API ", min_value=0.3)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.3)

    if st.button("Calculate Pb Petrosky"):
        x = 7.916*(10**-4)*(API**1.5410) - 4.561*(10**-5)*(T - 460)**1.3911
        Pb =((112.727*(Rs**.577421))/((Yg**.8439)*10**x))-1391.051
        st.success(f"Bubble point pressure Petrosky : {Pb:.4f}")
#..
def calculate_oil_formation_volume_factor_Standing():
    st.title("calculate oil formation volume Factor “Bo” using “Standing” Correlation")
    Rs = st.number_input("Enter the Gas solubility ", min_value=0.4)
    T = st.number_input("Enter the temperature (R)", min_value=0.4)
    Yo=st.number_input("Enter the oil specific gravity ", min_value=0.4)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.4)

    if st.button("Calculate Bo Standing"):
        Bo=.9759+.000120*(Rs*((Yg/Yo)**.5)+1.25*(T-460))**1.2
        st.success(f"oil formation volume standing : {Bo:.4f}")
#..
def calculate_oil_formation_volume_factor_Vasquez():
    st.title("calculate oil formation volume factor “Bo” using “Vasquez” Correlation")
    Rs = st.number_input("Enter the Gas solubility ", min_value=0.4)
    T = st.number_input("Enter the temperature (R)", min_value=0.4)
    API=st.number_input("Enter the API ", min_value=0.4)
    Tsep = st.number_input("Enter the temperature of the separator(R) ", min_value=0.4)
    Psep=st.number_input("Enter the pressure of the separator(psia) ", min_value=0.4)
    Yg=st.number_input("Enter the gas specific gravity ", min_value=0.4)

    if st.button("Calculate Bo Vasquez"):
        C1 = 4.677*10**-4 if API <= 30 else 4.670*10**-4
        C2 = 1.751*10**-5 if API <= 30 else 1.100*10**-5
        C3 = -1.811*10**-8 if API <= 30 else 1.337*10**-9
        Ygs = Yg*(1+(5.912*10**-5)*API*(Tsep-460)*math.log10(Psep/114.7))
        Bo = 1.0+C1*Rs+(T-520)*(API/Ygs)*(C2+C3*Rs)
        st.success(f"oil formation volume : {Bo:.4f}")
#..
def calculate_oil_formation_volume_factor_Glaso():
    st.title("calculate oil formation volume Factor “Bo” using “Glaso”  Correlation")
    Rs = st.number_input("Enter the Gas solubility ", min_value=0.4)
    T = st.number_input("Enter the temperature (R)", min_value=0.4)
    Yo=st.number_input("Enter the oil specific gravity ", min_value=0.4)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.4)

    if st.button("Calculate Bo Glaso"):
        Bob=Rs*(Yg/Yo)**.526+.968*(T-460)
        A=-6.58511+2.91329*math.log10(Bob)-.27683*(math.log10(Bob))**2
        Bo=1+10**A
        st.success(f"oil formation volume Glaso : {Bo:.4f}")
#..
def calculate_oil_formation_volume_factor_Marhoun():
    st.title("calculate oil formation volume Factor “Bo” using “Marhoun”  Correlation")
    Rs = st.number_input("Enter the Gas solubility ", min_value=0.4)
    T = st.number_input("Enter the temperature (R)", min_value=0.4)
    Yo=st.number_input("Enter the oil specific gravity ", min_value=0.4)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.4)

    if st.button("Calculate Bo Marhoun"):
        a=.742390
        b=.323294
        c=-1.202040
        F=(Rs**a)*(Yg**b)*(Yo**c)
        Bo=.497069+(.862963*10**-3)*T+(.182594*10**-2)*F+(.318009*10**-5)*F**2
        st.success(f"oil formation volume Marhoun : {Bo:.4f}")
#..

#..
def calculate_oil_formation_volume_factor_Petrosky():
    st.title("calculate oil formation volume Factor “Bo” using “Petrosky”  Correlation")
    Rs = st.number_input("Enter the Gas solubility ", min_value=0.4)
    T = st.number_input("Enter the temperature (R)", min_value=0.4)
    Yo=st.number_input("Enter the oil specific gravity ", min_value=0.4)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.4)

    if st.button("Calculate Bo Petrosky"):
        Bo=1.0113+(7.2046*10**-5)*((Rs**.3738)*(Yg**.2914/Yo**.6265)+.24626*(T-460)**.5371)**3.0936
        st.success(f"oil formation volume Petrosky : {Bo:.4f}")
#..
def calculate_oil_formation_volume_factor_materialbalance():
    st.title("calculate oil formation volume Factor “Bo” using “material balance”  equation")
    Rs = st.number_input("Enter the Gas solubility ", min_value=0.4)
    rho_o = st.number_input("Enter the density ", min_value=0.4)
    Yo=st.number_input("Enter the oil specific gravity ", min_value=0.4)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.4)

    if st.button("Calculate Bo material balance"):

        Bo=(62.4*Yo+.0136*Rs*Yg)/rho_o
        st.success(f"oil formation volume material balance : {Bo:.4f}")
#..
def calculate_crude_oil_density():
    st.title(" crude oil density “ρo” ")
    Rs = st.number_input("Enter the Gas solubility(scf/STB) ", min_value=0.4)
    Bo = st.number_input("Enter the oil formation volume factor ", min_value=0.4)
    Yo=st.number_input("Enter the oil specific gravity ", min_value=0.4)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.4)

    if st.button("Calculate ρo "):

        ρo=(62.4*Yo+.0136*Rs*Yg)/Bo
        st.success(f"crude oil density : {ρo:.4f}")
#..
def calculate_total_formation_volume_factor_standing():
    st.title("calculate total formation volume factor “Bt” using “standing”  correlation")
    Rs = st.number_input("Enter the Gas solubility(scf/STB) ", min_value=0.4)
    T = st.number_input("Enter the temperature (R) ", min_value=0.4)
    P=st.number_input("Enter the presssure (psia) ", min_value=0.4)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.4)
    Yo = st.number_input("Enter the oil specific gravity ", min_value=0.4)

    if st.button("Calculate Bt standing "):

        c=(2.9)*10**(-.00027*Rs)
        A_star=10**((math.log10(Rs*((T-460)**.5)*Yo**c)/Yg**.3)-(10.1-96.8/(6.604+math.log10(P))))
        Bt=10**(-5.223-(47.4/(-12.22+math.log10(A_star))))
        st.success(f"total formation volume factor(Bt) : {Bt:.4f}")
#..
def calculate_total_formation_volume_factor_Marhoun():
    st.title("calculate total formation volume factor “Bt” using “Marhoun”  correlation")
    Rs = st.number_input("Enter the Gas solubility(scf/STB) ", min_value=0.4)
    T = st.number_input("Enter the temperature (R) ", min_value=0.4)
    P=st.number_input("Enter the presssure (psia) ", min_value=0.4)
    Yg = st.number_input("Enter the gas specific gravity ", min_value=0.4)
    Yo = st.number_input("Enter the oil specific gravity ", min_value=0.4)

    if st.button("Calculate Bt Marhoun "):

        a=.644516
        b=-1.079340
        c=.724874
        d=2.006210
        e=-.761910
        F=(Rs**a)*(Yg**b)*(Yo**c)*(T**d)*(P**e)
        Bt=.314693+.106253*(10**-4)*F+.18883*(10**-10)*F**2
        st.success(f"total formation volume factor(Bt) : {Bt:.4f}")
#..
def calculate_crude_oil_viscosity_Dead_Beal():
    st.title("calculate crude DEAD oil viscosity “μod” using “Beal”  correlation")
    API = st.number_input("Enter the API ", min_value=0.4)
    T = st.number_input("Enter the temperature (R) ", min_value=0.4)


    if st.button("Calculate μod Beal "):

        a=10**(.43+(8.33/API))
        μod=(.32+(1.8*10**7)/(API**4.53))*(360/(T-260))**a

        st.success(f"crude oil viscosity(μod) : {μod:.4f}")
#..
def calculate_crude_oil_viscosity_Dead_Beggs():
    st.title("calculate crude DEAD oil viscosity “μod” using “Beggs”  correlation")
    API = st.number_input("Enter the API ", min_value=0.4)
    T = st.number_input("Enter the temperature (R) ", min_value=0.4)


    if st.button("Calculate μod Beggs "):

        Z=3.0324-.02023*API
        Y=10**Z
        X=Y*(T-460)**-1.163
        μod=(10**X)-1

        st.success(f"crude oil viscosity(μod) : {μod:.4f}")
#..
def calculate_crude_oil_viscosity_Dead_Glaso():
    st.title("calculate crude DEAD oil viscosity “μod” using “Glaso”  correlation")
    API = st.number_input("Enter the API ", min_value=0.4)
    T = st.number_input("Enter the temperature (R) ", min_value=0.4)


    if st.button("Calculate μod Glaso "):

        a=10.313*(math.log10(T-460))-36.447
        μod=(3.141*(10**10))*((T-460)**-3.444)*(math.log10(API))**a

        st.success(f"crude oil viscosity(μod) : {μod:.4f}")
#..
def calculate_crude_oil_viscosity_Saturated_Chew():
    st.title("calculate crude SATURATED oil viscosity “μob” using “Chew”  correlation")
    Rs = st.number_input("Enter the gas solubility ", min_value=0.4)
    μod = st.number_input("Enter the measured oil viscosity  ", min_value=0.4)


    if st.button("Calculate μob Chew "):

        e=3.74*(10**-3)*Rs
        d=1.1*(10**-3)*Rs
        c=8.62*(10**-5)*Rs
        b=(.68/10**c)+(.25/10**d)+(.062/10**e)
        a=Rs*(2.2*(10**-7)*Rs-7.4*(10**-4))
        μob=(10**a)*(μod)**b


        st.success(f"crude oil viscosity(μob) : {μob:.4f}")
#..
def calculate_crude_oil_viscosity_Saturated_Beggs():
    st.title("calculate crude SATURATED oil viscosity “μob” using  “Beggs” correlation")
    Rs = st.number_input("Enter the gas solubility ", min_value=0.4)
    μod = st.number_input("Enter the measured oil viscosity  ", min_value=0.4)


    if st.button("Calculate μob Beggs "):

        b=5.44*(Rs+150)**-.338
        a=10.715*(Rs+100)**-.515
        μob=a*(μod)**b
        st.success(f"crude oil viscosity(μob) : {μob:.4f}")
#..
def calculate_crude_oil_viscosity_underSaturated_Beal():
    st.title("calculate crude UNDER SATURATED oil viscosity “μo” using “Beal”  correlation")
    P = st.number_input("Enter the pressure ", min_value=0.4)
    Pb = st.number_input("Enter the pressure at bubble point  ", min_value=0.4)
    μob=st.number_input("Enter the Oil viscosity at the bubble-point pressure  ", min_value=0.0)


    if st.button("Calculate μo Beal "):

        μo=μob+.001*(P-Pb)*(.024*(μob**1.6)+.038*(μob**.56))
        st.success(f"crude oil viscosity(μo) : {μo:.4f}")
#..
def calculate_crude_oil_viscosity_underSaturated_Vasquez():
    st.title("calculate crude UNDER SATURATED oil viscosity “μo” using “Vasquez”  correlation")
    P = st.number_input("Enter the pressure ", min_value=0.4)
    Pb = st.number_input("Enter the pressure at bubble point  ", min_value=0.4)
    μob=st.number_input("Enter the Oil viscosity at the bubble-point pressure  ", min_value=0.0)


    if st.button("Calculate μo Vasquez "):

        a=-3.9*(10**-5)*P-5
        m=26*(P**1.187)*10**a
        μo=μob*(P/Pb)**m
        st.success(f"crude oil viscosity(μo) : {μo:.4f}")
#..
def calculate_Isothermal_Compressibility_Coefficient_of_Crude_Oil_vasquez ():
    st.title("calculate Isothermal Compressibility Coefficient of Crude Oil “co” using “Vasquez”  correlation")
    P = st.number_input("Enter the pressure (psia)", min_value=0.0)
    API=st.number_input("Enter the API gravity  ", min_value=0.0)
    T = st.number_input("Enter the Temperature  (R)", min_value=0.0)
    Rsb = st.number_input("Enter the gas solubility at bubble point pressure (scf/STB) ", min_value=0.0)
    Ygs = st.number_input("Enter the gas specific gravity at the reference separator pressure  ", min_value=0.0)


    if st.button("Calculate co Vasquez "):

        co=(-1433+5*17.2*Rsb*(T-460)-1180*Ygs+12.61*API)/(P*10**5)

        st.success(f"Isothermal Compressibility Coefficient of Crude Oil “co” : {co:.7f}")
#..
def calculate_Isothermal_Compressibility_Coefficient_of_Crude_Oil_Petrosky ():
    st.title("calculate Isothermal Compressibility Coefficient of Crude Oil “co” using “Petrosky”  correlation")
    P = st.number_input("Enter the pressure (psia)", min_value=0.0)
    API=st.number_input("Enter the API gravity  ", min_value=0.0)
    T = st.number_input("Enter the Temperature  (F)", min_value=0.0)
    Rsb = st.number_input("Enter the gas solubility at bubble point pressure (scf/STB) ", min_value=0.0)
    Yg = st.number_input("Enter the gas specific gravity   ", min_value=0.0)


    if st.button("Calculate co Petrosky "):

        co=(1.705*10**-7)*(Rsb**.6935)*(Yg**.1885)*(API**.3272)*(T**.6729)*P**-.5906

        st.success(f"Isothermal Compressibility Coefficient of Crude Oil “co” : {co:.7f}")
#..
def calculate_Isothermal_Compressibility_Coefficient_of_Crude_Oil_McCain ():
    st.title("calculate Isothermal Compressibility Coefficient of Crude Oil “co” using “McCain”  correlation")
    P = st.number_input("Enter the pressure (psia)", min_value=0.0)
    API=st.number_input("Enter the API gravity  ", min_value=0.0)
    T = st.number_input("Enter the Temperature  (R)", min_value=0.0)
    Rsb = st.number_input("Enter the gas solubility at bubble point pressure (scf/STB) ", min_value=0.0)
    Pb = st.number_input("Enter the bubble point pressure pressure (psia)  ", min_value=0.0)


    if st.button("Calculate co McCain "):

        A=-7.573-1.45*math.log(P)-.383*math.log(Pb)+1.402*math.log(T)+.256*math.log(API)+.449*math.log(Rsb)
        co=math.exp(A)

        st.success(f"Isothermal Compressibility Coefficient of Crude Oil “co” : {co:.7f}")
#..
def calculate_Isothermal_Compressibility_Coefficient_of_Crude_Oil ():
    st.title("calculate Isothermal Compressibility Coefficient of Crude Oil “co” using  correlation")
    P = st.number_input("Enter the pressure (psia)", min_value=0.0)
    Yo=st.number_input("Enter the oil specific gravity",format="%.3f", min_value=0.0)
    Yg = st.number_input("Enter the gas specific gravity",format="%.3f", min_value=0.0)
    T = st.number_input("Enter the Temperature  (R)", min_value=0.0)
    Rs = st.number_input("Enter the gas solubility , (scf/STB) ", min_value=0.0)
    Bo = st.number_input("oil formation volume factor at p, (bbl/STB)",format="%.3f", min_value=0.0)
    Bg = st.number_input("Enter the gas formation volume factor at pressure p (bbl/scf)",format="%.6f", min_value=0.0)


    if st.button("Calculate co  "):

        co=(-Rs/(Bo*(.83*P+21.75)))*(.00014*(Yg/Yo)**.5*((Rs*(Yg/Yo)**.5+1.25*(T-460))**.12)-Bg)

        st.success(f"Isothermal Compressibility Coefficient of Crude Oil “co” : {co:.7f}")
#..

#..

# Define the options for the select box
options = {
    'Calculate “γo” oil specific gravity': calculate_specific_gravity,
    'Calculate API gravity': calculate_API_gravity,
    'Calculate Rs using McCain Correlation':calculate_Gas_Solubility_McCain,
    'Calculate Rs using Standing Correlation':calculate_Gas_Solubility_standing,
    'Calculate Rs using Vasquez-Beggs correlation':calculate_gas_solubility_Vasquez,
    'Calculate Rs using Marhoun correlation':calculate_Gas_Solubility_Marhoun,
    'calculate Rs using Petrosky correlation':calculate_Gas_Solubility_Petrosky,
    'calculate Rs using Glaso correlation':calculate_Gas_Solubility_Glaso,
    'calculate Pb using standing correlation':calculate_Bubble_point_pressure_Standing,
    'calculate Pb using Vasquez correlation':calculate_Bubble_point_pressure_Vasquez,
    'calculate Pb using Glaso correlation':calculate_Bubble_point_pressure_Glaso,
    'calculate Pb using Marhoun correlation':calculate_Bubble_point_pressure_Marhoun,
    'calculate Pb using Petrosky correlation':calculate_Bubble_point_pressure_Petrosky,
    'calculate Bo using Standing correlation':calculate_oil_formation_volume_factor_Standing,
    'calculate Bo using Vasquez correlation':calculate_oil_formation_volume_factor_Vasquez,
    'calculate Bo using Glaso correlation':calculate_oil_formation_volume_factor_Glaso,
    'calculate Bo using Marhoun correlation':calculate_oil_formation_volume_factor_Marhoun,
    'calculate Bo using Petrosky correlation':calculate_oil_formation_volume_factor_Petrosky,
    'calculate Bo using materialbalance equation':calculate_oil_formation_volume_factor_materialbalance,
    'calculate ρo (crude oil density)':calculate_crude_oil_density,
    'calculate Bt using Standing correlation':calculate_total_formation_volume_factor_standing,
    'calculate Bt using Marhoun correlation':calculate_total_formation_volume_factor_Marhoun,
    'calculate μod  DEAD using Beal correlation':calculate_crude_oil_viscosity_Dead_Beal,
    'calculate μod  DEAD using Beggs correlation':calculate_crude_oil_viscosity_Dead_Beggs,
    'calculate μod  DEAD using Glaso correlation':calculate_crude_oil_viscosity_Dead_Glaso,
    'calculate μob  SATURATED using chew correlation':calculate_crude_oil_viscosity_Saturated_Chew,
    'calculate μob  SATURATED using Beggs correlation':calculate_crude_oil_viscosity_Saturated_Beggs,
    'calculate μo  UNDER SATURATED using Beal correlation':calculate_crude_oil_viscosity_underSaturated_Beal,
    'calculate μo  UNDER SATURATED using Vasquez correlation':calculate_crude_oil_viscosity_underSaturated_Vasquez,
    'calculate co  using vasquez correlation':calculate_Isothermal_Compressibility_Coefficient_of_Crude_Oil_vasquez,
    'calculate co  using petrosky correlation':calculate_Isothermal_Compressibility_Coefficient_of_Crude_Oil_Petrosky,
    'calculate co  using McCain correlation':calculate_Isothermal_Compressibility_Coefficient_of_Crude_Oil_McCain,
    'calculate co ':calculate_Isothermal_Compressibility_Coefficient_of_Crude_Oil


}

# Create the select box
selected_option = st.sidebar.selectbox('Select an option', list(options.keys()))

# Call the selected function
options[selected_option]()
#....

