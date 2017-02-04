/**
	main.cpp
	Usage: compile with C++11 or above. Run with ./a.out > image.ppm
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define f double
#define i int
#define RR return
#define O operator

#define P M_PI
#define M(x) (x>0)*x

// Vector class, 3 doubles, also used for storing color data in RGB format
struct vector {
	double x,y,z;
	vector(double a=0, double b=0, double c=0) {x=a; y=b; z=c;};
	vector operator+(vector o) {return {x+o.x, y+o.y, z+o.z};}  // Vector addition
	vector operator-(vector o) {return *this + o*-1.;}  // Vector subtraction
	vector operator*(double o) {return {o*x, o*y, o*z};}  // Multiplication by scalar
	double operator*(vector o) {return x*o.x + y*o.y + z*o.z;}  // Dot product
	vector operator!() {return *this*(1./sqrt(*this * *this));}  // Normalization
};
vector operator*(double a, vector o) {return o*a;}  // Left-multiplication by scalar

// Helper function, checks if a number is nearly-zero
bool z(double x) {
	return abs(x)<1e-4;
}

// Returns a random number between 0 and 1
double R() {
	return (double)rand()/RAND_MAX;
}

// Find nearest intersection, if any, between a ray with origin o and direction d and a sphere with center at the origin and unit radius. Returns t such that intersection = o + t*d. Assumes d is normalized
double B(vector o, vector d) {
	double t = -1*o*d;
	double l = o*o - t*t;
	if (l > 1)  // If l>1 no intersection, return 0
		return 0;
	double u = sqrt(1-l);
	double a = t - u;
	double b = t + u;

	// If a > 0 return a, else if b > 0 return b, else return 0 (fail, all intersections are behind)
	return a>0 ? a : (b>0 ? b : 0);
}

// Rotates a vector x with angle a around the x axis and angle b around the y axis
vector S(vector x, double a, double b) {
	return vector(
		        cos(b) * x.x + 0            +        sin(b) * x.z ,
		 sin(a)*sin(b) * x.x + cos(a) * x.y - sin(a)*cos(b) * x.z ,
	    -cos(a)*sin(b) * x.x + sin(a) * x.y + cos(a)*cos(b) * x.z 
	);
}

// Randomly rotates a vector a bit around its x and y axes
vector J(vector x, double e) {
	double t = (R()-.5)*P*e;
	double p = (R()-.5)*P*e;

	return S(x,t,p);
}

// The coordinates of the point light source
vector K = vector(3,2.5,2.5);

// Shoots a ray from point o with direction d (assumes d is normalized). N is the recursion counter, once it reaches 0 it doesn't recurse anymore (skips checking intersection vs sphere)
vector S(vector o, vector d, int N) {
	N--;  // Decrement the recursion counter
	vector p;

	// Compute first intersection with sphere
	double t = B(o, d);
	
	if (t>0 && N) {  // If there is a hit (and we haven't yet exhausted the recursion depth)
		vector p = o+t*d;
		vector n = !p;  // Normal
		
		// For the ball to appear slightly tilted, we want to rotate the frame of reference where we will make our calculations for the azimuthal and polar angles (to determine if the color of the point is yellow, blue or red inside the star) some angle about the x and y axes. This corresponds to the opposite rotation of the vector p.
		vector q = S(p,P/5,-P/3.5);
		
		// Here we calculate whether the point is inside the strip around the equator (blue) or not (yellow). Notice the calculations being done with respect to the rotated p
		vector C = abs(q.z)<cos(2*P/5.) ?
			vector(78, 191, 246) :
			vector(253, 252, 0);

		// Now we calculate whether the point is inside the star at the top. To do this we must know the polar angle (x), and we must do some geometric calculations. To derive these expressions write the straight lines for a 5-pointed star of minimum radius A and maximum radius B in polar coordinates around the origin
		double x = (atan2(q.y,q.x)+P)/2/P;
		double A = .16*P, B = .27*P;
		double a = A*cos(0.2*P)/(A*sin(0.2*P)-B);
		double b = tan(abs(fmod(x, 0.2)-0.1) * 2*P);
		double c = pow(a*B/(a-b),2)*(1+b*b);
		if (q.z > cos(c))
			C = vector(248,24,24);

		// Calculate reflected direction
		vector r = d-2*(d*n)*n;
		// Distance from point to the light
		vector L = K - p;
		
		// Lambertian shading, modified for the ilumination factor to be in the [.15,1] range
		C = C*(M(!L*!n)*.85+.15);

		// Reflections, recursively call S. Notice the bias on the origin vector to avoid surface acne, and the random jitter in the direction around the "true" direction for diffuse reflections
		vector G = S(p+1e-6*n, !J(r,0.05), N);
		
		// Specular highlights with the usual Phong shading formula
		vector H = pow(M(r*!L), 30)*255*vector(1,1,1);

		// Return a combination of the three effects
		return .82*C + .1*G + .25*H;

	} else {  // If it doesn't hit the sphere, it will hit a wall. Find out which wall
		vector W[] = {{-3,0,0},{0,9,0},{0,-2,0},{0,0,-1},{8,0,0},{0,0,8}};  // List of walls
		vector n;
		double u = 1e9;
		
		// Check point of intersection against each wall, if any
		for (auto w: W) {
			double t = d*w>1e-6 ? (w-o)*w/(d*w) : 1e9;
			// If an intersection is found, and it is closer than u, update u and p, as well as the normal
			if (t<u) {
				u = t;
				n = w;
				p = o+t*d;
			}
		}

		// Floor/wall colors
		vector a = {192,192,192};
		vector b = {184, 50, 50};
		
		// Vector from intersection point to light source, plus some jitter for soft shadows
		vector L = K + .15*vector(R(),R(),R()) - p;
		
		double I;  // Ilumination factor
		// Check if path to light is occluded by the sphere, if yes then it's a shadow
		double s = B(p, !L);
		if (s>0 && s*s<L*L) {
			I = 0.15;
		// Otherwise calculate illumination factor
		} else {
			I = M(-1.*!L*!n)*(0.5  // Lambertian shading
			   + 0.3*(1-tanh(16*acos(!L*!K) - 2.7*P)));  // Cone, for a "spotlight" effect (a Luxo lamp? ;))
		}

		// Subtract n from p to determine whether the wall is a side wall or floor/ceiling. Notice that whichever component of p-n is zero (or nearly-zero), it's going to be the non-zero component of the normal of the wall. If that is x or y, it's a side wall. If it's z, it's the floor or the ceiling.
		p = p-n;
		if (z(p.x)) {
			return (int(p.y*2+8)+int(p.z*2+8))&1 ? a*I : b*I;
		} else if (z(p.y)) {
			return (int(p.z*2+8)+int(p.x*2+8))&1 ? b*I : a*I;
		} else {
			vector r = !d-2*(!d*n)*n;
			return ((int(p.x*2+8)+int(p.y*2+8))&1 ? vector(28,28,28)*I : a*I) + .25*S(p+n, !J(r,0.08), N);
		}
	} 
}

int main() {
	// .ppm header (512x512 image with 2^8 channel depth)
	printf("P6 1920 1920 255 ");

	vector a = vector(6,.1,1),  // Position of the camera
		b = !vector(-1,0,0),  // Direction vector
		c = !vector(0,0,1),  // Up direction
		d = !vector(-.2,1,0);  // Right direction

	// Focal distance = distance from camera to the front of the ball, makes the ball be in focus
	double fd = 5.1;

	for (int y=1920; y--;) {  // Height
		for (int x=1920; x--;) {  // Width
			// Color starts at very dark grey, almost black
			vector p = {13,13,13};
			// Cast 64 rays
			int N = 64*4;
			for (int r=N; r--;) {
				// Random delta to sum to starting position, to simulate nonzero aperture (for depth of field blur)
				vector t = .1*(c*(R()-.5) + d*(R()-.5));
				// Cast a ray from starting position (pinhole camera position + random delta), in the correct direction for the current pixel (with random deltas for stochastic sampling), and subtract the t/fd factor to keep focal plane in focus. Accumulate onto p
				p = p + 1./N * S(a+t, !(b + 0.00125/2*(d*(480*2-x+R())-c*(480*2-y+R()))-1./fd*t ), 8);
			}
			// Clamp values to 255
			if (p.x >= 255) p.x = 255;
			if (p.y >= 255) p.y = 255;
			if (p.z >= 255) p.z = 255;
			// Print pixel values
			printf("%c%c%c", int(p.x), int(p.y), int(p.z));
		}
	}
}