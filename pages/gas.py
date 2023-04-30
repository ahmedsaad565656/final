import streamlit as st
import math
def calculate_gas_specific_gravity():
    st.title("Gas Specific Gravity “γg”")


    # Ask the user for the inputs
    Ma = st.number_input("Enter the apparent molecular weight of the gas (Ma) ", min_value=0.0)

    # Calculate the gas specific gravity and display the result
    if st.button("Calculate γg"):
        γg = Ma/28.96
        st.success(f"gas specific gravity: {γg:.4f}")
#..
def calculate_gas_specific_volume():
    st.title("Gas Specific volume “v”")


    # Ask the user for the inputs
    ρg = st.number_input("Enter the gas density (ρg), (lb/ft3)  ", min_value=0.0)

    # Calculate the gas specific gravity and display the result
    if st.button("Calculate v"):
        v=1/ρg
        st.success(f"gas specific volume: {v:.4f} ft3/lb" )
#..
def calculate_gas_viscosity_standing():
    st.title(" CORRECTED VISCOSITY OF NATURAL GASES“μg” using “standing” correlation ")


    # Ask the user for the inputs
    Yg = st.number_input("Enter the gas specific gravity (Yg)  ", min_value=0.0)
    T = st.number_input("Enter the reservoir Temperature (T), (in R)  ", min_value=0.0)
    yco2 = st.number_input("Enter the mole fraction of CO2 (yco2)  ", min_value=0.0)
    yN2 = st.number_input("Enter the mole fraction of N2 (yco2)  ", min_value=0.0)
    yH2S = st.number_input("Enter the mole fraction of H2S (yco2)  ", min_value=0.0)

    # Calculate the gas specific gravity and display the result
    if st.button("Calculate μg standing "):
        μCo2=yco2*((9.08*10**-3)*math.log10(Yg)+6.24*10**-3)
        μN2 = yN2 * ((8.48 * 10 ** -3) * math.log10(Yg) + 9.59 * 10 ** -3)
        μH2s = yH2S * ((8.49 * 10 ** -3) * math.log10(Yg) + 3.73 * 10 ** -3)
        μuncorrected=(1.709*10**-5-2.062*10**-6*Yg)*(T-460)+8.118*10**-3-6.15*10**-3*math.log10(Yg)
        μ= μuncorrected+μCo2+μN2+μH2s

        st.success(f"standing viscosity: {μ:.4f} cp" )
#..
def calculate_gas_viscosity_lee():
    st.title("VISCOSITY OF NATURAL GASES“μg” using “lee” correlation ")


    # Ask the user for the inputs
    Ma = st.number_input("Enter the apparent molecular weight of gas mixture (Ma)  ", min_value=0.0)
    T = st.number_input("Enter the reservoir pressure (T), (in R)  ", min_value=0.0)
    ρg = st.number_input("Enter the gas density (ρg), (in lb/ft3)  ", min_value=0.0)

    # Calculate the gas specific gravity and display the result
    if st.button("Calculate μg lee"):
        X=3.5+(986/T)+.01*Ma
        Y=2.4-.2*X
        K=((9.4+.02*Ma)*T**1.5)/(209+19*Ma+T)
        μg=(K*10**-4)*math.exp(X*((ρg/62.4)**Y))

        st.success(f"lee viscosity: {μg:.4f} cp" )
#..
def calculate_gas_viscosity_Dempsey():
    st.title("VISCOSITY OF NATURAL GASES“μg” using “Dempsey” correlation ")


    # Ask the user for the inputs
    Pr = st.number_input("Enter the pseudo-reduced pressure of the gas mixture (Pr), in psia  ", min_value=0.0)
    Tr = st.number_input("Enter the pseudo-reduced temperature of the gas mixture (Tr), (in R)  ", min_value=0.0)
    μ1 = st.number_input("Enter the  viscosity of the gas at atmospheric pressure and reservoir temperature,(μ1) (in cp)  ", min_value=0.0)


    # Calculate  and display the result
    if st.button("Calculate μg Dempsey"):
        a0=-2.46211820
        a1=2.970547414
        a2=-2.86264054*10**-1
        a3=8.05420522*10**-3
        a4=2.80860949
        a5=-3.49803305
        a6=3.60373020*10**-1
        a7=-1.044324*10**-2
        a8=-7.93385648*10**-1
        a9=1.39643306
        a10=-1.49144925*10**-1
        a11=4.41015512*10**-3
        a12=8.39387178*10**-2
        a13=-1.86408848*10**-1
        a14=2.03367881*10**-2
        a15=-6.09579263*10**-4

        X=a0+a1*Pr+a2*Pr**2+a3*Pr**3+Tr*(a4+a5*Pr+a6*Pr**2+a7*Pr**3)+Tr**2*(a8+a9*Pr+a10*Pr**2+a11*Pr**3)+Tr**3*(a12+a13*Pr+a14*Pr**2+a15*Pr**3)-math.log(Tr,math.exp(1))
        μg=μ1*math.exp(X)



        st.success(f"Dempsey viscosity: {μg:.4f} cp" )
#..

#..
def calc_gas_formation_volume_factor_Z():
    st.title("calculate gas formation volume factor “Bg” in terms of z factor")


    # Ask the user for the inputs
    units = st.selectbox("Select the units for the gas formation volume factor", ("ft^3/scf", "bbl/scf"))
    p = st.number_input("Enter the gas pressure (in psia)", min_value=0.0)
    T = st.number_input("Enter the gas temperature (in °R)", min_value=0.0)
    Z = st.number_input("Enter the gas compressibility factor", min_value=0.0)


    # Calculate the gas specific gravity and display the result
    if st.button("Calculate Bg"):
        if units == "ft^3/scf":
            # Calculate gas formation volume factor using the formula Bg = 0.02827*Z*T/p (in ft^3/scf)
            Bg = 0.02827 * Z * T / p
        elif units == "bbl/scf":
            # Calculate gas formation volume factor using the formula Bg = 0.005035*Z*T/p (in bbl/scf)
            Bg = 0.005035 * Z * T / p
        st.success(f"gas formation volume factor: {Bg:.4f} {units}")
#..
def calc_gas_formation_volume_factor_density():
    st.title("calculate gas formation volume factor “Bg” in terms of density")


    # Ask the user for the inputs
    units = st.selectbox("Select the units for the gas formation volume factor", ("ft^3/scf", "bbl/scf"))
    ρg = st.number_input("Enter the gas density (ρg),(in lb/ft3)", min_value=0.0)
    Ma = st.number_input("Enter the gas  apparent molecular weight (Ma)", min_value=0.0)


    # Calculate the gas specific gravity and display the result
    if st.button("Calculate Bg"):
        if units == "ft^3/scf":
            # Calculate gas formation volume factor using the formula Bg = 0.02827*Z*T/p (in ft^3/scf)
            Bg = 0.002635*(Ma/ρg)
        elif units == "bbl/scf":
            # Calculate gas formation volume factor using the formula Bg = 0.005035*Z*T/p (in bbl/scf)
            Bg = 0.000469*Ma/ρg
        st.success(f"gas formation volume factor : {Bg:.4f} {units}")
#..
def calc_gas_expansion_volume_factor_Z():
    st.title("calculate gas expansion volume factor “Eg” in terms of z factor")


    # Ask the user for the inputs
    units = st.selectbox("Select the units for the gas formation volume factor", ("ft^3/scf", "bbl/scf"))
    p = st.number_input("Enter the gas pressure (in psia)", min_value=0.0)
    T = st.number_input("Enter the gas temperature (in °R)", min_value=0.0)
    Z = st.number_input("Enter the gas compressibility factor", min_value=0.0)


    # Calculate the gas specific gravity and display the result
    if st.button("Calculate Eg"):
        if units == "ft^3/scf":
            # Calculate gas formation volume factor using the formula Bg = 0.02827*Z*T/p (in ft^3/scf)
            Eg = 1/(0.02827 * Z * T / p)
        elif units == "bbl/scf":
            # Calculate gas formation volume factor using the formula Bg = 0.005035*Z*T/p (in bbl/scf)
            Eg = 1/(0.005035 * Z * T / p)
        st.success(f"gas expansion volume factor: {Eg:.4f} {units}")
#..
def calc_gas_expansion_volume_factor_density():
    st.title("calculate gas expansion volume factor “Eg” in terms of density")


    # Ask the user for the inputs
    units = st.selectbox("Select the units for the gas formation volume factor", ("ft^3/scf", "bbl/scf"))
    ρg = st.number_input("Enter the gas density (ρg),(in lb/ft3)", min_value=0.0)
    Ma = st.number_input("Enter the gas  apparent molecular weight (Ma)", min_value=0.0)


    # Calculate the gas specific gravity and display the result
    if st.button("Calculate Eg"):
        if units == "ft^3/scf":
            # Calculate gas formation volume factor using the formula Bg = 0.02827*Z*T/p (in ft^3/scf)
            Eg = 1/(0.002635*(Ma/ρg))
        elif units == "bbl/scf":
            # Calculate gas formation volume factor using the formula Bg = 0.005035*Z*T/p (in bbl/scf)
            Eg = 1/(0.000469*Ma/ρg)
        st.success(f"gas expansion volume factor : {Eg:.4f} {units}")
#..
def calc_compressibility_of_natural_gases():
    st.title("calculate compressibility of natural gases “Cg” ")


    # Ask the user for the inputs
    type = st.selectbox("Select the type of the gas", ("Ideal", "Real"))
    if type == "Ideal":
       P=st.number_input("Enter the pressure (p),(in psia)", min_value=0.0)
    elif type == "Real":
        P = st.number_input("Enter the pressure (p),(in psia)", min_value=0.0)
        Pr = st.number_input("Enter the Pseudo reduced pressure (Pr)", min_value=0.0)
        Tr = st.number_input("Enter the Pseudo reduced temperature(Tr)", min_value=0.0)


    # Calculate the gas specific gravity and display the result
    if st.button("Calculate cg"):
        if type == "Ideal":
            # Calculate gas formation volume factor using the formula Bg = 0.02827*Z*T/p (in ft^3/scf)
            Cg = 1/P
        elif type == "Real":
            # Calculate gas formation volume factor using the formula Bg = 0.005035*Z*T/p (in bbl/scf)
            z = 1.008505 + 0.04623 * (Pr / Tr) + (0.862707 * Pr ** 1.368627) / (10 ** (0.636778 * Tr)) - (
                    2.324825 * Pr) / (10 ** (0.649787 * Tr))
            d = 0.036662 / Tr + 1.12921 * Pr ** 0.60015 / 10 ** (0.783224 * Tr) - 2.950267 * Pr / 10 ** (
                    0.789306 * Tr)
            Cg = 1 / P - d / z
        st.success(f"compressibility of natural gas : {Cg:.4f} psi^-1")
#..
def calc_z_factor_Tarek():
    st.title("Calculate gas compressibility factor z using Tarek Ahmed correlation “Z” ")

    Pr = st.number_input("Enter the Pseudo reduced pressure (Pr)", min_value=0.0)
    Tr = st.number_input("Enter the Pseudo reduced temperature(Tr)", min_value=0.0)

    if st.button("Calculate z"):
        z = 1.008505 + 0.04623 * (Pr / Tr) + (0.862707 * Pr ** 1.368627) / (10 ** (0.636778 * Tr)) - (2.324825 * Pr) / (
                    10 ** (0.649787 * Tr))

        st.success(f"Compressibility of natural gas: {z:.4f}")
#..
def calc_z_factor_Tarek():
    st.title("Calculate gas compressibility factor “Z” using Tarek Ahmed correlation  ")

    Pr = st.number_input("Enter the Pseudo reduced pressure (Pr)", min_value=0.0)
    Tr = st.number_input("Enter the Pseudo reduced temperature(Tr)", min_value=0.0)

    if st.button("Calculate z"):
        z = 1.008505 + 0.04623 * (Pr / Tr) + (0.862707 * Pr ** 1.368627) / (10 ** (0.636778 * Tr)) - (2.324825 * Pr) / (
                    10 ** (0.649787 * Tr))

        st.success(f"gas compressibility factor: {z:.4f}")
#..
def calc_z_factor_Hall():
    st.title("Calculate gas compressibility factor “Z” using Hall correlation  ")

    Pr = st.number_input("Enter the Pseudo reduced pressure (Pr)", min_value=0.0)
    T = st.number_input("Enter the  temperature(T),(in R)", min_value=0.0)
    Tc = st.number_input("Enter the Pseudo critical temperature (Pc)", min_value=0.0)

    if st.button("Calculate z"):
        t=Tc/T
        Y=.0125*Pr*t*math.exp(-1.2*(1-t)**2)
        z=(.06125*Pr*t/Y)*math.exp(-1.2*(1-t)**2)

        st.success(f"gas compressibility factor: {z:.4f}")
#..
def calc_z_factor_Hall():
    st.title("Calculate gas compressibility factor “Z” using Hall correlation  ")

    Pr = st.number_input("Enter the Pseudo reduced pressure (Pr)", min_value=0.0)
    T = st.number_input("Enter the  temperature(T),(in R)", min_value=0.0)
    Tc = st.number_input("Enter the Pseudo critical temperature (Pc),(in R)", min_value=0.0)

    if st.button("Calculate z"):
        t=Tc/T
        Y=.0125*Pr*t*math.exp(-1.2*(1-t)**2)
        z=(.06125*Pr*t/Y)*math.exp(-1.2*(1-t)**2)

        st.success(f"gas compressibility factor: {z:.4f}")
#..
def calc_z_factor_DAK():
    st.title("Calculate gas compressibility factor “Z” using DAK correlation  ")

    Pr = st.number_input("Enter the Pseudo reduced pressure (Pr)", min_value=0.0)
    ρr = st.number_input("Enter the  reduced gas density ρr", min_value=0.0)
    Tr = st.number_input("Enter the Pseudo reduced temperature (Pr)", min_value=0.0)

    if st.button("Calculate z"):

        z=(.27*Pr)/(Tr*ρr)

        st.success(f"gas compressibility factor: {z:.4f}")
#..
def calc_z_factor_Papay ():
    st.title("Calculate gas compressibility factor “Z” using Papay  correlation  ")

    Pr = st.number_input("Enter the Pseudo reduced pressure (Pr)", min_value=0.0)

    Tr = st.number_input("Enter the Pseudo reduced temperature (Pr)", min_value=0.0)

    if st.button("Calculate z"):

        z=1-(3.53*Pr)/(10**(.9813*Tr))+((.274*Pr**2)/10**(.8157*Tr))

        st.success(f"gas compressibility factor: {z:.4f}")
#..
def calc_z_factor_Gopal ():
    st.title("Calculate gas compressibility factor “Z” using Gopal  correlation  ")

    Pr = st.number_input("Enter the Pseudo reduced pressure (Pr)", min_value=0.0)

    Tr = st.number_input("Enter the Pseudo reduced temperature (Pr)", min_value=0.0)

    if st.button("Calculate z"):

        z=Pr*(.711+3.66*Tr)**-1.4667-1.6371*(.319*Tr+.522)+2.071

        st.success(f"gas compressibility factor: {z:.4f}")
#..
def calc_pseudo_critical_pressure ():
    st.title("Calculate pseudo critical pressure  “Ppc” ")
    gas_type = st.selectbox("choose the type of gas system", ("Natural Gas Systems", "Gas-Condensate Systems"))
    Yg = st.number_input("Enter the gas specific gravity (Yg)", min_value=0.0)
    Pc = None

    if st.button("Calculate Pc"):
        if gas_type == "Natural Gas Systems":
            Pc = 667 + 15 * Yg - 37.5 * Yg ** 2
        elif gas_type == "Gas-Condensate Systems":
            Pc = 706 + 51.7 * Yg - 11.1 * Yg ** 2

    if Pc is not None:
        st.success(f"pseudo critical pressure: {Pc:.4f} R")
#..
def calc_pseudo_critical_temperature ():
    st.title("Calculate pseudo critical temperature  “Tpc” ")
    gas_type = st.selectbox("choose the type of gas system", ("Natural Gas Systems", "Gas-Condensate Systems"))
    Yg = st.number_input("Enter the gas specific gravity (Yg)", min_value=0.0)
    Tc = None

    if st.button("Calculate Tc"):
        if gas_type == "Natural Gas Systems":
            Tc = 168+325*Yg-12.5*Yg**2
        elif gas_type == "Gas-Condensate Systems":
            Tc = 187+330*Yg-71.5*Yg**2

    if Tc is not None:
        st.success(f"pseudo critical temperature: {Tc:.4f} R")
#..
def calc_corrected_pseudo_critical_temperature_wichert ():
    st.title(" calculate corrected pseudo critical temperature using wichert correlation  ")



    Tc = st.number_input("Enter the Pseudo reduced temperature (Tc),(in R)", min_value=0.0)
    YH2S = st.number_input("Enter the mole fraction of H2S in the gas mixture (YH2S)", min_value=0.0)
    YCO2 = st.number_input("Enter the Enter the mole fraction of CO2 in the gas mixture (YCO2)", min_value=0.0)

    if st.button("Calculate Tc corrected"):

        A=YH2S+YCO2
        E=120*(A**.9-A**1.6)+15*(YH2S**.5-YH2S**4)
        T_corrected=Tc-E


        st.success(f"corrected pseudo critical temperature: {T_corrected:.4f}")
#..
def calc_corrected_pseudo_critical_pressure_wichert ():
    st.title(" calculate corrected pseudo critical pressure using wichert correlation  ")

    Pc = st.number_input("Enter the Pseudo critical pressure (Pc),(in psia)", min_value=0.0)
    Tc = st.number_input("Enter the Pseudo critical temperature (Tc),(in R)", min_value=0.0)
    Tc_corrected = st.number_input("Enter the  corrected Pseudo critical temperature (Tc_corrected),(in R)", min_value=0.0)
    YH2S = st.number_input("Enter the mole fraction of H2S in the gas mixture (YH2S)", min_value=0.0)
    YCO2 = st.number_input("Enter the Enter the mole fraction of CO2 in the gas mixture (YCO2)", min_value=0.0)

    if st.button("Calculate Pc corrected"):

        A=YH2S+YCO2
        E=120*(A**.9-A**1.6)+15*(YH2S**.5-YH2S**4)
        p_corrected=(Pc*Tc_corrected)/(Tc+YH2S*(1-YH2S)*E)


        st.success(f"corrected pseudo critical pressure: {p_corrected:.4f}")
#..
def calc_corrected_pseudo_critical_temperature_carr ():
    st.title(" calculate corrected pseudo critical temperature using carr correlation  ")


    Tr = st.number_input("Enter the  corrected Pseudo reduced temperature (Tr_corrected),(in R)", min_value=0.0)
    YH2S = st.number_input("Enter the mole fraction of H2S in the gas mixture (YH2S)", min_value=0.0)
    YCO2 = st.number_input("Enter the Enter the mole fraction of CO2 in the gas mixture (YCO2)", min_value=0.0)
    YN2 = st.number_input("Enter the Enter the mole fraction of N2 in the gas mixture (YN2)", min_value=0.0)

    if st.button("Calculate Tc corrected"):

        Tc_corrected=Tr-80*YCO2+130*YH2S-250*YN2



        st.success(f"corrected pseudo critical temperature: {Tc_corrected:.4f} R")
#..
def calc_corrected_pseudo_critical_pressure_carr ():
    st.title(" calculate corrected pseudo critical pressure using carr correlation  ")


    Pc = st.number_input("Enter the  corrected Pseudo critical pressure (Pc_corrected),(in psia)", min_value=0.0)
    YH2S = st.number_input("Enter the mole fraction of H2S in the gas mixture (YH2S)", min_value=0.0)
    YCO2 = st.number_input("Enter the Enter the mole fraction of CO2 in the gas mixture (YCO2)", min_value=0.0)
    YN2 = st.number_input("Enter the Enter the mole fraction of N2 in the gas mixture (YN2)", min_value=0.0)

    if st.button("Calculate Pc corrected"):

        Pc_corrected=Pc+440*YCO2+600*YH2S-170*YN2



        st.success(f"corrected pseudo critical pressure: {Pc_corrected:.4f} psia")
#..
def calc_apparent_molecular_weight ():
    st.title(" calculate the apparent molecular weight  “Ma”")

    num_components = st.number_input("Enter the number of components in the gas mixture", min_value=1, step=1)
    yi = []
    Mi = []
    for i in range(num_components):
        yi_i = st.number_input(f"Enter the mole fraction (yi) for component {i + 1}", min_value=0.0, max_value=1.0,
                               step=0.01)
        Mi_i = st.number_input(f"Enter the molecular weight (Mi) for component {i + 1}", min_value=0.0)
        yi.append(yi_i)
        Mi.append(Mi_i)

    if st.button("Calculate Ma"):

        Ma = sum([yi[i]*Mi[i] for i in range(len(yi))])



        st.success(f"Apparent Molecular weight: {Ma:.4f} ")
#..
options = {
    'Calculate  “Ma”': calc_apparent_molecular_weight,
    'Calculate “γg” gas specific gravity': calculate_gas_specific_gravity,
    'Calculate “v” gas specific volume': calculate_gas_specific_volume,
    'Calculate “μg” using standing correlation': calculate_gas_viscosity_standing,
    'Calculate “μg” using lee correlation': calculate_gas_viscosity_lee,
    'Calculate “μg” using  Dempsey correlation': calculate_gas_viscosity_Dempsey,
    'Calculate “Bg” in terms of Z ': calc_gas_formation_volume_factor_Z,
    'Calculate “Bg” in terms of density ': calc_gas_formation_volume_factor_density,
    'Calculate “Eg” in terms of Z ': calc_gas_expansion_volume_factor_Z,
    'Calculate “Eg” in terms of density ': calc_gas_expansion_volume_factor_density,
    'Calculate “Cg”  ': calc_compressibility_of_natural_gases,
    'Calculate “z” using tarek correlation  ': calc_z_factor_Tarek,
    'Calculate “z” using Hall correlation  ': calc_z_factor_Hall,
    'Calculate “z” using DAK correlation  ': calc_z_factor_DAK,
    'Calculate “z” using Papay correlation  ': calc_z_factor_Papay,
    'Calculate “z” using Gopal correlation  ': calc_z_factor_Gopal,
    'Calculate “Pc”   ': calc_pseudo_critical_pressure,
    'Calculate “Tc”   ': calc_pseudo_critical_temperature,
    'Calculate “Tc corr ” wichert   ': calc_corrected_pseudo_critical_temperature_wichert,
    'Calculate “Pc corr”  wichert ': calc_corrected_pseudo_critical_pressure_wichert,
    'Calculate “Tc corr”  carr': calc_corrected_pseudo_critical_temperature_carr,
    'Calculate “Pc corr”  carr': calc_corrected_pseudo_critical_pressure_carr,


}

# Create the select box
selected_option = st.sidebar.selectbox('Select an option', list(options.keys()))

# Call the selected function
options[selected_option]()
#..

