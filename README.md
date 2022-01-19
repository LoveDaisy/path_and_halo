# Path and halo

It is a simple try. It is quite straight to trace rays from sun to crystal, and trace every refraction and
reflection in crystal, and finally we get the halo image.

And what if we do this in a backword way? We collect all possible ray path through a crystal, and when we want to
check and render halo, we pick a right ray path from the pre-computed bank and check its intensity,
then we set it as the intensity of the point in the image.

That is the motivation of this project.

## Quick start

I use a thirdparty library to deal with EALPix grids, see [acknowledgement](#acknowledgement).
So you must put path of MEALPix to matlab system path in order to use them.

`test_packages.m` makes some basic tests.

`main.m` set up a simple configuration for a crystal and sun, and run the calculation.

## Acknowledgement

I use a matlab version HEALPix algorithm, called MEALPix. It is forked from https://sourceforge.net/projects/mealpix/
