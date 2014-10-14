// Functions for the planets

// Copyright Ole Nielsen 2002-2004
// Please read copyright notice in astrotools2.html source

// Formulae and elements from Paul Schlyter's article "Computing planetary positions" available at 
// http://hem.passagen.se/pausch/comp/ppcomp.html

MERCURY = 0; VENUS = 1; EARTH = 2; MARS = 3; JUPITER = 4; SATURN = 5; 
URANUS = 6; NEPTUNE = 7; SUN = 9; MOON = 10; COMET = 15; USER = 20;

// Planet diameters at 1 AU in arcsec (km for Moon)
var ndiam = [6.72, 16.68, 1, 9.36, 196.88, 165.46, 70.04, 67.0, 1, 1919.3, 716900000.0];

// The planet object

function planet(name,num,N,i,w,a,e,M) {
	this.name=name;
	this.num=num;
	this.N=N; 	// longitude of ascending node
	this.i=i;		// inclination
	this.w=w;	// argument of perihelion
	this.a=a;	// semimajor axis
	this.e=e;		// eccentricity
	this.M=M;	// mean anomaly
}

// elements from Paul Schlyter
var planets=new Array();
planets[0]=new planet("Mercury",0,
   new Array(48.3313, 3.24587E-5),
   new Array(7.0047, 5.00E-8),
   new Array(29.1241, 1.01444E-5),
   new Array(0.387098, 0),
   new Array(0.205635, 5.59E-10),
   new Array(168.6562, 4.0923344368));

planets[1]=new planet("Venus  ",1,
   new Array(76.6799, 2.46590E-5),
   new Array(3.3946, 2.75E-8),
   new Array(54.8910, 1.38374E-5),
   new Array(0.723330, 0),
   new Array(0.006773, -1.302E-9),
   new Array(48.0052, 1.6021302244));

planets[2]=new planet("Earth  ",2,
   new Array(0,0),
   new Array(0,0),
   new Array(0,0),
   new Array(0.0, 0.0),
   new Array(0.0, 0.0),
   new Array(0,0));

planets[3]=new planet("Mars   ",3,
   new Array(49.5574, 2.11081E-5),
   new Array(1.8497, -1.78E-8),
   new Array(286.5016, 2.92961E-5),
   new Array(1.523688, 0),
   new Array(0.093405, 2.516E-9),
   new Array(18.6021, 0.5240207766));

planets[4]=new planet("Jupiter",4,
   new Array(100.4542, 2.76854E-5),
   new Array(1.3030, -1.557E-7),
   new Array(273.8777, 1.64505E-5),
   new Array(5.20256, 0),
   new Array(0.048498, 4.469E-9),
   new Array(19.8950, 0.0830853001));

planets[5]=new planet("Saturn ",5,
	new Array(113.6634, 2.38980E-5),
	new Array(2.4886, -1.081E-7),
	new Array(339.3939, 2.97661E-5),
	new Array(9.55475, 0),
	new Array(0.055546, -9.499E-9),
	new Array(316.9670, 0.0334442282));

planets[6]=new planet("Uranus ",6,
	new Array(74.0005, 1.3978E-5),
	new Array(0.7733, 1.9E-8),
	new Array(96.6612, 3.0565E-5),
	new Array(19.18171, -1.55E-8),
	new Array(0.047318, 7.45E-9),
	new Array(142.5905, 0.011725806));

planets[7]=new planet("Neptune",7,
	new Array(131.7806, 3.0173E-5),
	new Array(1.7700, -2.55E-7),
	new Array(272.8461, -6.027E-6),
	new Array(30.05826, 3.313E-8),
	new Array(0.008606, 2.15E-9),
	new Array(260.2471, 0.005995147));

// Body holds current data of planet, Sun or Moon), method .update(jday,obs)
function Body(name,number,colour,colleft,colright) {
	this.name=name;
	this.number=number;
	this.colour=colour;
	this.colleft=colleft;
	this.colright=colright;
	this.alt=0;
	this.az=0;
	this.dec=0;
	this.ra=0;
	this.H=0;
	this.eclon=0;
	this.eclat=0;
	this.illum=1;
	this.r=1;	// heliocentric distance
	this.dist=1;	//geocentric distance
	this.mag=-1.0;
	this.elong=0;
	this.pa=0;	// position angle (elongation)
	this.update=updatePosition;	
	this.elongupdate=updateElong;
}

bodies = new Array();
bodies[0] = new Body("Mercury",MERCURY,0,24,25);
bodies[1] = new Body("Venus  ",VENUS,1,24,25);
bodies[2] = new Body("Earth  ",2,3,24,25);
bodies[3] = new Body("Mars   ",MARS,3,24,25);
bodies[4] = new Body("Jupiter",JUPITER,4,24,25);
bodies[5] = new Body("Saturn ",SATURN,5,24,25);
bodies[6] = new Body("Uranus ",URANUS,6,24,25);
bodies[7] = new Body("Neptune",NEPTUNE,7,24,25);
bodies[8] = new Body("",8,0,0,0);
bodies[9] = new Body("Sun    ",SUN,9,26,27);
bodies[10] = new Body("Moon   ",MOON,10,28,29);
bodies[COMET] = new Body("Comet  ",COMET,2,24,25);
bodies[20] = new Body("User object",USER,2,24,25);

function updatePosition(jday,obs) {
// update body-object with current positions
// elongation NOT calculated! (use updateElong for that)
	this.p=this.number;
	if (this.p==USER) {		// fixed/user object
		var altaz=radec2aa(this.ra, this.dec, jday, obs);
		this.alt=altaz[0];
		this.az=altaz[1];
		this.H=altaz[2];
		return;
	}
	var dat=PlanetAlt(this.p,jday,obs);
	this.alt=dat[0];
	this.az=dat[1];
	this.H=dat[2];
	this.ra=dat[3];
	this.dec=dat[4] - (dat[4]>180.0 ? 360 : 0);
	this.eclon=rev(dat[5]);
	this.eclat=dat[6];
	this.r=dat[8];
	this.dist=dat[9];
	this.illum=dat[7];
	this.mag=dat[10];
}

function updateElong(jday,obs) {
	// calculate elongation for object
	if (this.number==SUN) return;
	bodies[SUN].update(jday,obs); 
	var ra2=bodies[SUN].ra; var dec2=bodies[SUN].dec;
	this.update(jday,obs);
	var dat = separation(this.ra, ra2, this.dec, dec2);
	this.elong = dat[0];
	this.pa = dat[1];
}

// heliocentric xyz for planet p
// this is not from Meeus' book, but from Paul Schlyter http://hem.passagen.se/pausch/comp/ppcomp.html
// account for pertuberations of Jupiter, Saturn, Uranus (Uranus and Neptune mutual pertubs are included in elements)
// returns heliocentric x, y, z, distance, longitude and latitude of object
function helios(p,jday) {
	var d = jday-2451543.5;
	var w = p.w[0] + p.w[1]*d;		// argument of perihelion
	var e = p.e[0] + p.e[1]*d; 
	var a = p.a[0] + p.a[1]*d;
	var i = p.i[0] + p.i[1]*d;
	var N = p.N[0] + p.N[1]*d;
	var M = rev( p.M[0] + p.M[1]*d ); 	// mean anomaly
	var E0 = M + RAD2DEG*e*sind(M) * ( 1.0+e*cosd(M) );
	var E1 = E0 - ( E0-RAD2DEG*e*sind(E0)-M ) / ( 1.0-e*cosd(E0) );
	while (Math.abs(E0-E1) > 0.0005) {
		E0 = E1;
		E1 = E0 - ( E0 - RAD2DEG*e*sind(E0)-M ) / ( 1.0-e*cosd(E0) );
	}
	var xv = a*(cosd(E1) - e);
	var yv = a*Math.sqrt(1.0 - e*e) * sind(E1);
	var v = rev(atan2d( yv, xv ));		// true anomaly
	var r = Math.sqrt( xv*xv + yv*yv );	// distance
	var xh = r * ( cosd(N)*cosd(v+w) - sind(N)*sind(v+w)*cosd(i) );
	var yh = r * ( sind(N)*cosd(v+w) + cosd(N)*sind(v+w)*cosd(i) );
	var zh = r * ( sind(v+w)*sind(i) );
	var lonecl = atan2d(yh, xh);
	var latecl = atan2d(zh, Math.sqrt(xh*xh + yh*yh + zh*zh));
	if (p.num==JUPITER) {		// Jupiter pertuberations by Saturn
		var Ms = rev(planets[SATURN].M[0] + planets[SATURN].M[1]*d);
		lonecl += (-0.332)*sind(2*M-5*Ms-67.6) - 0.056*sind(2*M-2*Ms+21) + 0.042*sind(3*M-5*Ms+21) -
				0.036*sind(M-2*Ms) + 0.022*cosd(M-Ms) + 0.023*sind(2*M-3*Ms+52) - 0.016*sind(M-5*Ms-69);
		xh=r*cosd(lonecl)*cosd(latecl);		// recalc xh, yh
		yh=r*sind(lonecl)*cosd(latecl);
	}
	if (p.num==SATURN) {		// Saturn pertuberations
		var Mj = rev(planets[JUPITER].M[0] + planets[JUPITER].M[1]*d);
		lonecl += 0.812*sind(2*Mj-5*M-67.6) - 0.229*cosd(2*Mj-4*M-2) + 0.119*sind(Mj-2*M-3) + 
				0.046*sind(2*Mj-6*M-69) + 0.014*sind(Mj-3*M+32);
		latecl += -0.020*cosd(2*Mj-4*M-2) + 0.018*sind(2*Mj-6*M-49);
		xh = r*cosd(lonecl)*cosd(latecl);		// recalc xh, yh, zh
		yh = r*sind(lonecl)*cosd(latecl);
    	zh = r*sind(latecl);
	}
	if (p.num==URANUS) {		// Uranus pertuberations
		var Mj = rev(planets[JUPITER].M[0] + planets[JUPITER].M[1]*d);
		var Ms = rev(planets[SATURN].M[0] + planets[SATURN].M[1]*d);
		lonecl += 0.040*sind(Ms-2*M+6) + 0.035*sind(Ms-3*M+33) - 0.015*sind(Mj-M+20);
		xh=r*cosd(lonecl)*cosd(latecl);		// recalc xh, yh
		yh=r*sind(lonecl)*cosd(latecl);
	}
	return new Array(xh,yh,zh,r,lonecl,latecl);
}	// helios()


function radecr(obj,sun,jday,obs) {
// radecr returns ra, dec and earth distance
// obj and sun comprise Heliocentric Ecliptic Rectangular Coordinates
// (note Sun coords are really Earth heliocentric coordinates with reverse signs)
		// Equatorial geocentric co-ordinates
		var xg=obj[0]+sun[0];
		var yg=obj[1]+sun[1];
		var zg=obj[2];
		// Obliquity of Ecliptic (exponent corrected, was E-9!)
		var obl = 23.4393 - 3.563E-7 * (jday-2451543.5);
		// Convert to eq. co-ordinates
		var x1=xg;
		var y1=yg*cosd(obl) - zg*sind(obl);
		var z1=yg*sind(obl) + zg*cosd(obl);
		// RA and dec (33.2)
		var ra=rev(atan2d(y1, x1));
		var dec=atan2d(z1, Math.sqrt(x1*x1+y1*y1));
		var dist=Math.sqrt(x1*x1+y1*y1+z1*z1);
		return new Array(ra,dec,dist);
} 


function radec2aa(ra,dec,jday,obs) {
// Convert ra/dec to alt/az, also return hour angle, azimut = 0 when north
// DOES NOT correct for parallax!
// TH0=Greenwich sid. time (eq. 12.4), H=hour angle (chapter 13)
		var TH0 = 280.46061837 + 360.98564736629*(jday-2451545.0);
		var H = rev(TH0-obs.longitude-ra);
		var alt = asind( sind(obs.latitude)*sind(dec) + cosd(obs.latitude)*cosd(dec)*cosd(H) );
		var az = atan2d( sind(H), (cosd(H)*sind(obs.latitude) - tand(dec)*cosd(obs.latitude)) );
		return new Array (alt, rev(az+180.0), H);
}


function separation(ra1, ra2, dec1, dec2) {		
// ra, dec may also be long, lat, but PA is relative to the chosen coordinate system
	var d = acosd(sind(dec1)*sind(dec2) + cosd(dec1)*cosd(dec2)*cosd(ra1-ra2));		// (Meeus 17.1)
	if (d < 0.1) d = Math.sqrt(sqr( rev2(ra1-ra2)*cosd((dec1+dec2)/2) ) + sqr(dec1-dec2));	// (17.2)
	var pa = atan2d(sind(ra1-ra2),cosd(dec2)*tand(dec1)-sind(dec2)*cosd(ra1-ra2));		// angle
	return new Array(d, rev(pa));
}	// end separation()


function PlanetAlt(p,jday,obs) {
// Alt/Az, hour angle, ra/dec, ecliptic long. and lat, illuminated fraction, dist(Sun), dist(Earth), brightness of planet p
		if (p==SUN) return SunAlt(jday,obs);
		if (p==MOON) return MoonPos(jday,obs);
		if (p==COMET) return CometAlt(jday,obs);
		var sun_xyz = sunxyz(jday);
		var planet_xyz = helios(planets[p],jday);

		var dx = planet_xyz[0]+sun_xyz[0];
		var dy = planet_xyz[1]+sun_xyz[1];
		var dz = planet_xyz[2]+sun_xyz[2];
		var lon = rev( atan2d(dy, dx) );
		var lat = atan2d(dz, Math.sqrt(dx*dx+dy*dy));

		var radec = radecr(planet_xyz, sun_xyz, jday, obs);
		var ra = radec[0]; 
		var dec = radec[1];
		var altaz = radec2aa(ra, dec, jday, obs);

		var dist = radec[2];
		var R = sun_xyz[3];	// Sun-Earth distance
		var r = planet_xyz[3];	// heliocentric distance
		var k = ((r+dist)*(r+dist)-R*R) / (4*r*dist);		// illuminated fraction (41.2) 
		// brightness calc according to Meeus p. 285-86 using Astronomical Almanac expressions
		var absbr = new Array(-0.42, -4.40, 0, -1.52, -9.40, -8.88, -7.19, -6.87);
		var i = acosd( (r*r+dist*dist-R*R) / (2*r*dist) );	// phase angle
		var mag = absbr[p] + 5 * log10(r*dist);	// common for all planets
		switch(p) {
		case MERCURY:
			mag += i*(0.0380 + i*(-0.000273 + i*0.000002)); break;
		case VENUS:
			mag += i*(0.0009 + i*(0.000239 - i*0.00000065)); break;
		case MARS:
			mag += i*0.016; break;
		case JUPITER:
			mag += i*0.005; break;
		case SATURN:	// (Ring system needs special treatment, see Meeus Ch. 45)
			var T = (jday-2451545.0)/36525;		// (22.1)
			var incl = 28.075216 - 0.012998*T + 0.000004*T*T;	// (45.1)
			var omega = 169.508470 + 1.394681*T + 0.000412*T*T;	// (45.1)
			var B = asind(sind(incl)*cosd(lat)*sind(lon-omega) - cosd(incl)*sind(lat));
			var l = planet_xyz[4];	// heliocentric longitude of Saturn
			var b = planet_xyz[5];	// heliocentric latitude (do not confuse with 'b' in step 6, page 319)
			// correction for Sun's aberration skipped
			var U1 = atan2d(sind(incl)*sind(b)+cosd(incl)*cosd(b)*sind(l-omega), cosd(b)*cosd(l-omega));
			var U2 = atan2d(sind(incl)*sind(lat)+cosd(incl)*cosd(lat)*sind(lon-omega), cosd(lat)*cosd(lon-omega));
			var dU = Math.abs(U1 - U2);
			mag += 0.044*dU - 2.60*sind(Math.abs(B)) + 1.25*sind(B)*sind(B);
			break;
		}
		return new Array (altaz[0], altaz[1], altaz[2], ra, dec, lon, lat, k, r, dist, mag); 
}

