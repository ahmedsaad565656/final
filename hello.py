import streamlit as st

st.set_page_config(
    page_title="Hello",
    page_icon="👋",
)

st.write("# Welcome to Reservoir Properties calculator! 👋")

st.sidebar.success("Select a category above.")

st.markdown(
    """
    A reservoir properties calculator is an application that can be used to estimate various physical properties of hydrocarbon reservoirs such as oil, gas, and water. This type of application typically takes input data such as reservoir pressure, temperature, and fluid composition, and then uses mathematical models and empirical correlations to calculate properties such as reservoir fluid densities, viscosities, and phase behavior.
    **👈 Select a category from the sidebar** to see the properties of each category!
    ### check the references
    -  Reservoir engineering handbook By Tarek Ahmed fifth edition.[Here](https://www.academia.edu/40740457/Reservoir_Engineering_Handbook_Ahmed_Tarek_5th_edition)
    - A Study on Gas Compressibility Factor for Gas-Condensate Systems . [Here](https://jpme.journals.ekb.eg/article_38796_d2bc02bd3e0f1f9d7e732f288f2f5f4d.pdf)
    - State of the Art – Natural Gases Viscosity under Reservoir Conditions .[Here](https://onepetro.org/SPESATS/proceedings-abstract/05TSSA/All-05TSSA/SPE-106326-MS/188071)
    ### contacts
    - I hope that this reservoir properties calculator app is helpful and provides accurate and reliable estimates of reservoir properties for oil, gas, and water. If you encounter any issues or have suggestions for improvements, please do not hesitate to contact me. Your feedback is important to me and will help to ensure that this app is as useful and effective as possible. Thank you for using this app and for your support!
    
    - you can find me here [LinkedIn](linkedin.com/in/ahmed-saad-b5b211203/),[Facebook](https://www.facebook.com/profile.php?id=100054258981827)
"""
)