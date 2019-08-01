Back to [Fan](https://fanwangecon.github.io)'s
[Dynamic Assets Repository](https://fanwangecon.github.io/CodeDynaAsset/) Table of Content.

# The ABZ model with Borrowing

- *a*: asset state/choice
- *b*: borrowing, this version of the code allows for borrowing, with and without default
- *z*: shock

# Descriptions

See descriptions for **Section 2** on the [Dynamic Assets Repository](https://fanwangecon.github.io/CodeDynaAsset/) main page.

# Effects of Changing Savings and Borrowing Interest rates

Among other simulations, we simulate:

1. *abz* the effects **borrowing** parameters on distributional outcomes
    * cross test, **no default**, [adjust borrow rate, bounds, etc](https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_nbc_cross.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_abz/test/ff_az_ds_vecsv/test_borr/fsi_abz_ds_vecsv_nbc_cross.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_nbc_cross.html)
    * cross test, **default**, [adjust borrow rate, bounds, etc](https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default_cross.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_abz/test/ff_az_ds_vecsv/test_borr/fsi_abz_ds_vecsv_default_cross.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default_cross.html)    
2. *abz* the effects **min inc and save r** on outcomes
    * cross test, **no default**, [min income and savings r](https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_price/html/fsi_abz_ds_vecsv_price_nbc_cross.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_abz/test/ff_az_ds_vecsv/test_price/fsi_abz_ds_vecsv_price_nbc_cross.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_price/html/fsi_abz_ds_vecsv_price_nbc_cross.html)
    * cross test, **default**, [min income and savings r](https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_price/html/fsi_abz_ds_vecsv_price_default_cross.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_abz/test/ff_az_ds_vecsv/test_price/fsi_abz_ds_vecsv_price_default_cross.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_price/html/fsi_abz_ds_vecsv_price_default_cross.html)    


Note that interest rate changes have price and substitution effects, but additionally, like other parameters, they will also impact the distribution.

## What is the effect of increasing savings interest rate

As one would expect, people are "better off" when the savings interest rate increases, in terms of their consumption and savings.

Without default, as savings interest rates increase, aggregate savings go up, consumption go up. Standard savings increases. Standard deviation of consumption could decrease than increase. With default, mean consumptions are slightly higher in generally at all savings interest rate levels.

With default, seems to have slightly lower consumption variance, and also slightly higher consumption mean. Perhaps makes sense, default allows for wider range of borrowing, so better consumption smoothing?
