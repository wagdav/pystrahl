from nose.tools import assert_raises, assert_equal

import strahl


def test_sanity_check():
    p = strahl.default_params()

    strahl.sanity_check(p)  # this should go fine

    p.pop('shot')  # remove a parameter
    assert_raises(AssertionError, strahl.sanity_check, p)


from strahl.datablocks import plasma_background, geometry, main_parameter_file
def test_datablocks():
    p = strahl.default_params()

    plasma_background(p)
    geometry(p)
    main_parameter_file(p)


from strahl.create_input import input_directory_names, input_file_names
def test_create_input():
    o = input_directory_names('xxx')
    assert_equal(o['nete'], 'xxx/nete')
    assert_equal(o['param_files'], 'xxx/param_files')
    assert_equal(o['results'], 'xxx/results')

    o = input_file_names('xxx', strahl.default_params())
    assert_equal(o['geometry'], 'xxx/nete/grid_99999.0')
    assert_equal(o['background'], 'xxx/nete/pp99999.0')
    assert_equal(o['main'], 'xxx/param_files/main_99999.0')
    assert_equal(o['influx'], 'xxx/nete/Arflx99999.dat')
