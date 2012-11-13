from nose.tools import assert_raises

import strahl


def test_sanity_check():
    p = strahl.default_params()

    strahl.sanity_check(p)  # this should go fine

    p.pop('shot')  # remove a parameter
    assert_raises(AssertionError, strahl.sanity_check, p)
