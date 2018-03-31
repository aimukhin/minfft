# fourier

## Implementation details
Here is a short explanation of how real transforms are reduced to a
half-length complex transform.

Everywhere below,

![](exp.svg)

### Real DFT
Results of the real DFT can be recovered as

![](realdft-recover.svg)

from the results of the complex DFT of length N/2:

![](dft-for-realdft.svg)

Inverse real DFT routine simply plays backwards the steps taken by the
forward transform routine.

### Type-2 real symmetric transforms
DCT-2 and DST-2 of length N can be easily reduced to a real DFT of the
same length:

![](type2-impl.svg)

### Type-3 real symmetric transforms
Since the type-3 transforms are inverses (up to a multiple) of their
type-2 counterparts, they are computed by working backwards the steps
taken by the type-2 routines.

### Type-4 real symmetric transforms
These transforms can be written as

![](type4-impl.svg)

Consider

![](g-def.svg)

and notice symmetry

![](g-symmetry.svg)

Therefore, we will determine all G if we find those at the even
positions. Luckily, there's a simple way to compute them:

![](g-even.svg)

This way we reduce a type-4 transform to a half-length complex DFT.
