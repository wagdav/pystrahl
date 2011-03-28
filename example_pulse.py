import strahl

params = strahl.defaultParams()

t, flx = strahl.rectangular_pulse(length=5e-3, max_value=2.5e23)

params['impurity.influx'] = (t, flx)
params['numerical.time.dt'] = 1e-4
params['numerical.time.final'] = 0.5

strahl.create_input(params)
