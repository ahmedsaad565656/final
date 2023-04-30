import streamlit as st
import math

def calculate_water_formation_volume_factor_free():
    st.title("calculate water formation volume factor“Bw” (gas free water)")
    P = st.number_input("Enter the pressure (psia)", min_value=0.0)
    T = st.number_input("Enter the temperature (R)", min_value=0.0)

    if st.button("Calculate Bw free"):
        A1 = .9947 + (5.8 * 10 ** -6) * (T - 460) + (1.02 * 10 ** -6) * (T - 460) ** 2
        A2 = (-4.228 * 10 ** -6) + (1.8376 * 10 ** -8) * (T - 460) + (-6.77 * 10 ** -11) * (T - 460) ** 2
        A3 = (1.3 * 10 ** -10) + (-1.3855 * 10 ** -12) * (T - 460) + (4.285 * 10 ** -15) * (T - 460) ** 2
        Bw = A1 + A2 * P + A3 * (P ** 2)
        st.success(f"water formation volume factor Bw : {Bw:.4f}")

def calculate_water_formation_volume_factor_saturated():
    st.title("calculate water formation volume factor “Bw”(gas saturated water)")
    P = st.number_input("Enter the pressure (psia)", min_value=0.0)
    T = st.number_input("Enter the temperature (R)", min_value=0.0)

    if st.button("Calculate Bw saturated"):
        A1 = .9911 + (6.35 * 10 ** -5) * (T - 460) + (8.5 * 10 ** -7) * (T - 460) ** 2
        A2 = (-1.093 * 10 ** -6) + (-3.497 * 10 ** -9) * (T - 460) + (-4.57 * 10 ** -12) * (T - 460) ** 2
        A3 = (-5* 10 ** -11) + (6.429 * 10 ** -13) * (T - 460) + (-1.43 * 10 ** -15) * (T - 460) ** 2
        Bw = A1 + A2 * P + A3 * (P ** 2)
        st.success(f"water formation volume factor Bw : {Bw:.4f}")
#..
def calculate_water_Isothermal_Compressibility():
    st.title("calculate water Isothermal Compressibility “Cw”")
    P = st.number_input("Enter the pressure (psia) ", min_value=0.0)
    T = st.number_input("Enter the temperature (F)", min_value=0.0)

    if st.button("Calculate Cw "):
        c1=3.8546-.000134*P
        c2=-.01052+(4.77*10**-7)*P
        c3=(3.9267*10**-5)-(8.8*10**-10)*P
        Cw=(c1+c2*T+c3*T**2)*10**-6
        st.success(f"water Isothermal Compressibility  : {Cw:.10f} psi^-1")
#..
def calculate_Gas_Solubility_in_Water():
    st.title("calculate Gas Solubility in Water “Rsw”")
    P = st.number_input("Enter the pressure (psia) ", min_value=0.0)
    T = st.number_input("Enter the temperature (F)", min_value=0.0)

    if st.button("Calculate Rsw "):
        A=2.12+3.45*(10**-3)*T-3.59*(10**-5)*T**2
        B=.0107-5.26*(10**-5)*T+1.48*(10**-7)*T**2
        C=8.75*(10**-7)+3.9*(10**-9)*T-1.02*(10**-11)*T**2
        Rsw=A+B*P+C*P**2
        st.success(f"Gas Solubility in Water  : {Rsw:.4f} ")
#..
def calculate_water_viscosity_Brill():
    st.title("calculate water viscosity “μw” using “Brill” method")

    T = st.number_input("Enter the temperature (R)", min_value=0.0)

    if st.button("Calculate μw Brill "):
        μw=math.exp(1.003-(1.479*10**-2)*(T-460)+(10**-5)*(T-460)**2)
        st.success(f"water viscosity  : {μw:.4f} cp")
#..
def calculate_water_viscosity_Meehan():
    st.title("calculate water viscosity “μw” using “Meehan” method")
    P = st.number_input("Enter the pressure (psia)", min_value=0.0)
    T = st.number_input("Enter the temperature (R)", min_value=0.0)
    Ws = st.number_input("Enter the weight percent of salt in brine ", min_value=0.0)

    if st.button("Calculate μw Meehan"):
        D=1.12166-.0263951*Ws+(6.79461*(10**-4)*Ws**2)+(5.47119*(10**-5)*Ws**3)-(1.55586*(10**-6)*Ws**4)
        μwt=(109.574-8.40564*Ws+(.313314*(Ws**2))+(8.72213*(10**-3)*Ws**3))*(T-460)**-D
        μw=μwt*(.9994+(1.0295*10**-5)*P+(3.1062*10**-9)*P**2)
        st.success(f"water viscosity  : {μwt:.4f} cp")
#..
def calculate_water_density():
    st.title("calculate water density “ρw” ")

    Ws = st.number_input("Enter the weight percent of salt in brine ", min_value=0.0)

    if st.button("Calculate ρw"):
        ρw=62.368+.438603*Ws+(1.60074*10**-9)*Ws**2
        st.success(f"water density  : {ρw:.4f} ")
# Define the options for the select box
options = {

    'Calculate  Bw (gas free water) ': calculate_water_formation_volume_factor_free,
    'Calculate  Bw (gas saturated water) ': calculate_water_formation_volume_factor_saturated,
    'Calculate  Cw  ': calculate_water_Isothermal_Compressibility,
    'Calculate  Rsw  ': calculate_Gas_Solubility_in_Water,
    'Calculate μw using Brill correlation ': calculate_water_viscosity_Brill,
    'Calculate μw using Meehan correlation': calculate_water_viscosity_Meehan,
    'Calculate ρw': calculate_water_density

}

# Create the select box
selected_option = st.sidebar.selectbox('Select an option', list(options.keys()))

# Call the selected function
options[selected_option]()
#.................................

