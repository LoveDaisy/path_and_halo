# Path and halo

It is a simple try.

As we know, it is quite straight to trace rays from sun to crystal, and trace every refraction and
reflection in crystal, and finally we get the halo image.

And what if we do this in a vice versa way?
Say, if we are able to collect/compute and store all possible ray path through a crystal in a bank,
then when we want to query a intensity at a given direction, we simply filter out those raypaths
that make rays along this direction, and sum them up to get the intensity.

That is the motivation of this project.

I originally write in MATLAB, and I am working on re-write in C++. Check folder `matlab` and `cpp`
for detail.
