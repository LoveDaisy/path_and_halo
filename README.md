# Path and halo

It is a simple try. I find it intuitive to trace rays from sun to crystal, and trace every refraction and reflection in crystal, and finally we get the halo image.

And what if we do this in a backword way? We compute all possible ray path through a crystal, and when we want to check and render halo, we pick a right ray path from the pre-computed bank and check its intensity, then we set it as the intensity of the point in the image.

## Acknowledgement

I use a matlab version HEALPix algorithm, called MEALPix. It is forked from https://sourceforge.net/projects/mealpix/
