import strahl

rc = strahl.defaults.defaultParams.copy()

t, flx = strahl.rectangular_pulse(5e-3, 2.5e23)

rc['impurity.influx.time'] = t
rc['impurity.influx.flux'] = flx
rc['numerical.time.dt'] = 1e-4
rc['numerical.time.final'] = 0.5

strahl.create_input(rc)
