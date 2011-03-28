import strahl

rc = strahl.defaultParams()

t, flx = strahl.rectangular_pulse(lenght=5e-3, max_value=2.5e23)

rc['impurity.influx'] = (t, flx)
rc['numerical.time.dt'] = 1e-4
rc['numerical.time.final'] = 0.5

strahl.create_input(rc)
