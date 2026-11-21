# COJ en05 hybrid order solutions (AVA B60)

en05 (COJ FLOYDS camera until early 2019) has no sky flats in the LCO archive. Each file here is a hybrid order
solution built by `banzai_floyds_offline en05-orders` (`banzai_floyds/en05.py`) from one en05 lamp flat:
order 1 traced on the lamp flat, order 2 the en12 sky-flat shape of the same slit width (2022-03-13 1.6",
2022-03-22 2.0" / 6.0") shifted onto the lamp light. Solutions from the same year agree to < 0.1 px rms.

| file | slit | order-2 shift |
|---|---|---|
| coj2m002-en05-20140902-0006-f91 | 1.6 | +2.42 |
| coj2m002-en05-20140902-0009-f91 | 2.0 | +1.40 |
| coj2m002-en05-20150705-0012-f91 | 2.0 | +3.60 |
| coj2m002-en05-20160616-0003-f91 | 2.0 | +12.68 |
| coj2m002-en05-20170602-0003-f91 | 6.0 | +8.90 |
| coj2m002-en05-20170612-0007-f91 | 2.0 | +11.70 |
| coj2m002-en05-20180607-0016-f91 | 2.0 | +8.02 |

They have no window in `skyflats.dat`, so banzai uses the nearest one in time with the frame's slit width.
`setup_db` registers them; for an existing calibration database run `banzai_floyds_offline register-en05`.
