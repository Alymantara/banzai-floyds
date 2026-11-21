"""AVA B60: COJ en05 (no sky flats) gets en12's order/LSF configuration and hybrid order solutions."""
import numpy as np
from astropy.table import Table

from banzai_floyds import en05


def _coeffs(c0_1, c0_2, dom2=(615.0, 1920.0)):
    return Table({'c0': [c0_1, c0_2], 'c1': [1.0, 2.0], 'c2': [0.1, 0.2], 'c3': [0.0, 0.0], 'c4': [0.0, 0.0],
                  'order': [1, 2], 'domainmin': [55.0, dom2[0]], 'domainmax': [1600.0, dom2[1]], 'height': [93, 93]})


def test_hybrid_coeffs_keeps_order1_and_shifts_reference_order2():
    own, reference = _coeffs(120.0, 999.0, (1250.0, 1920.0)), _coeffs(124.0, 400.0)
    new = en05.hybrid_coeffs(own, reference, 2.5)
    assert list(new[0]) == list(own[0])                                    # order 1: the lamp-flat trace
    assert new['c0'][1] == 402.5 and new['c1'][1] == 2.0                   # order 2: reference shape + shift
    assert (new['domainmin'][1], new['domainmax'][1]) == (615.0, 1920.0)   # and its full domain


def test_register_en05_copies_en12_configuration(tmp_path):
    from banzai.dbs import get_session
    from banzai_floyds import offline
    from banzai_floyds.dbs import LSFParams, OrderHeight, OrderLocation
    address = f"sqlite:///{tmp_path / 'floyds.db'}"
    offline.setup_db(address)                                              # calls register_en05
    assert en05.register_en05(address, solutions=[]) == 0                  # again: idempotent
    e12, e05 = en05._instrument(address, 'en12'), en05._instrument(address, 'en05')
    assert e05 is not None and e05.site == 'coj' and e05.id != e12.id
    with get_session(db_address=address) as session:
        for table in (OrderHeight, OrderLocation, LSFParams):
            n12 = session.query(table).filter(table.instrument_id == e12.id).count()
            n05 = session.query(table).filter(table.instrument_id == e05.id).count()
            assert n12 > 0 and n05 == n12, table.__tablename__


def test_fit_heights_keeps_orders_on_the_detector():
    """B60: en05's blue order (centre y ~ 470 of 512) gets a smaller height; an order with room keeps 93."""
    t = _coeffs(120.0, 470.0)
    t['c1'] = [10.0, 0.0]; t['c2'] = [0.0, 0.0]
    new = en05.fit_heights(t, (512, 2048))
    assert new['height'][0] == 93                                          # order 1: y 110..130, plenty of room
    assert new['height'][1] == 2 * (511 - 470 - 2)                         # order 2: 78, so 470 + 39 <= 509
