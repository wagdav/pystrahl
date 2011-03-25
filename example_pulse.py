import strahl
import defaults

bg = defaults.background.copy()
main = defaults.main.copy()

bg['influx'] = strahl.rectangular_pulse(5e-3, 2.5e23)
main['dt'] = 1e-4
main['tfinal'] = 0.5

strahl.create_input(defaults.geometry, bg, main)
