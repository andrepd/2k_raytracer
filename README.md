"Pixar Ball" Raytracer in 2kB
=============================

Inspired by the famous "[minimal raytracer challenge][1]" by Paul S. Heckbert, I wrote this program as a challenge to myself to see how much I could cram into a raytracer that was at most 2kB of source code.

`main.cpp` weighs exactly 2000 bytes. `pretty.cpp` is the same code, only properly formatted and commented.

Features
--------
* Textured yellow Pixar ball with blue stripe around the equator and red star on top, placed at an angle
* Glossy effect on ball with blurred reflections
* Ambient spherical light + cone spotlight
* Phong highlights
* Textured floor and walls
* Reflective floor
* Soft shadows
* Depth of field blur

Output
------
You should get an output like this (click on image to see a larger version):

[![Sample output](http://img.pixady.com/2017/02/129628_hdthumb2.png)](http://img.pixady.com/2017/02/686553_hd.png)

Usage
-----
Compile with `g++ -O3 -std=c++11 main.cpp`. File is output to stdout in the PPM file format. Run with `./a.out > img.ppm`. It should complete in less than a minute on a modern machine. Open with a compatible viewer.

License
------
This work is in the public domain.

![Public domain](https://i.creativecommons.org/p/mark/1.0/88x31.png)

[1]: https://books.google.ca/books?id=CCqzMm_-WucC&pg=PA375&lpg=PA375&dq=%22A+Minimal+Ray+Tracer%22&source=bl&ots=msmz42NHge&sig=rYEdHlY0zC2Sk_aPaZhzjMhyfj8&hl=en&sa=X&ei=2NQ8Utb2I-ae2gWHn4GYAg#v=onepage&q=%22A%20Minimal%20Ray%20Tracer%22&f=false
