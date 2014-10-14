// This Node Module repurposes code created by Ole Nielsen and Peter Hayes                                           
// Copyright Ole Nielsen 2002-2004
// Copyright Peter Hayes 1999-2001
// Published under GNU License
// Jeff Gladchun 2014


exports.returnAllPlanets = function(obs) {
	// debugger;
	// print data for all planets for ONE day. Positions are for local noon.
	var obscopy=new Object();
	for (var i in obs) {
		obscopy[i] = obs[i];
	}
	// obscopy.hours = 0;	// set to local midnight
	// obscopy.minutes = 0;

	var ndiam = [6.72, 16.68, 1, 9.36, 196.88, 165.46, 70.04, 67.0, 1, 1919.3, 716900000.0];

	var objects=[9,10,0,1,3,4,5,6,7];

	var result = {
		observer: obs,
		coords: []
	};

	for (var n in objects) {
		var objData = {};
		var obj = objects[n]
		var jday = jd(obscopy);

		bodies[obj].update(jday, obs);
		bodies[obj].elongupdate(jday, obs);
		
		var pa = bodies[obj].pa;

		objData.name = bodies[obj].name.trim();
		objData.rightAscension = hmsstring(bodies[obj].ra/15);
		objData.declination = anglestring(bodies[obj].dec,false,true);
		objData.elong = (fixnum(bodies[obj].elong,6,1) + "° " + (pa >= 180 ? " W" : " E")).trim();
		objData.dist = bodies[obj].dist;

		result.coords.push(objData);
	}
	return result;
};

function nextDate(obs1,dstep,origday) {
	// update date and time, origday is day of month at start (some months may not allow this day)
	if (dstep < 0) { 	// dstep is in months (integer!)
		dstep = -dstep;
		if (dstep >= 12) {
			obs1.year += Math.floor(dstep/12); 
		}
		else {
			obs1.month += dstep;
			if (obs1.month > 12) {
				obs1.year++;
				obs1.month -= 12;
			}
		}
		month_length[1]=leapyear(obs1.year)?29:28;	// check for leapyear
		obs1.day = ( origday>month_length[obs1.month-1] ? month_length[obs1.month-1] : origday);
	}
	else {	// dstep is in days (max 31)
		var m = Math.round(1440*(dstep-Math.floor(dstep)));
		obs1.minutes += m - 60*(Math.floor(m/60));
		obs1.hours += Math.floor(m/60);
		obs1.day += Math.floor(dstep);
		if (obs1.minutes > 59) {
			obs1.minutes -= 60;
			obs1.hours++;
		}
		if (obs1.hours > 23) {
			obs1.hours -= 24;
			obs1.day++;
		}
		month_length[1]=leapyear(obs1.year)?29:28;	// check for leapyear
		while (obs1.day > month_length[obs1.month-1]) {
			obs1.day -= month_length[obs1.month-1];
			obs1.month++; 
			if (obs1.month == 13) {
				obs1.year++;
				obs1.month = 1;
			}
		}
	}
}		// end nextDate()


function dfrac2tstr(t) {
	// returns time string from fraction of day (0 <= t < 1). If t < 0 return '--:--'
	if (t < 0 || t >= 1) return "--:--";
	t1 = Math.round(1440*t);	// round to nearest minute
	var hours = Math.floor(t1/60);
	var minutes = t1 - 60*hours;
	return hmstring2(hours, minutes,0);
}


function pheader(doc,obj,obs,title,descrip,line1,line2) {
	// common code for page header
	var str = head1 + "Javascript: " + title + head2;
	str += "<p><A HREF=\"javascript:window.close()\">close window</A></p>\n";
	str += "<h2>" + title + "</h2><p><b>" + descrip + "</b></p>";
	if (obj>=0 && obj<100) str += "<h3>Object: " + bodies[obj].name + "</h3>";
	if (obj==100) str += "<h3>Object: All planets</h3>";
	str += "<p>Observer location: "+sitename();
	str += " (UT "+hmstring(-obs.tz/60.0,true)+")</p>\n";
	var line3="";
	for (var i=0; i<line2.length; i++) line3 += "-";
	str += "<pre>" + line1 + "\n" + line2 + "\n" + line3 + "\n";
	doc.write(str);
}		//	end pheader()


function pbottom(doc,pwin,line2) {
	// finish the page
	var line3="";
	for (var i=0; i<line2.length; i++) line3 += "-";
	var str = line3 + "</pre>\n";
	str += "<p><A HREF=\"javascript:window.close()\">close window</A></p>\n";
	str += "</CENTER></BODY></HTML>";
	doc.write(str);
	doc.close();
	pwin.focus();
}		// end pbottom()


function doTwilightAltAz(obs,obj,dspan,dstep,sunalt) {
	// Altitude of object at specified Sun depression (twilight visibility)
	var obscopy=new Object(); var obsmax=new Object();
	for (var i in obs) {
		obscopy[i] = obs[i]; obsmax[i] = obs[i];
	}
	obscopy.hours = 12;	// set to local noon
	obscopy.minutes = 0;
	if (dstep > 0 && dstep < 1.0) dstep = 1.0;	// time step must be at least one day

	var pwin=window.open("","twilight","menubar,scrollbars,resizable");
	var doc=pwin.document;
	var descrip="at the times when the Sun is " + (-sunalt) + "&deg; below the horizon";
	var line1="                     Morning Twilight                    Evening Twilight  ";
	var line2="   Date          Time     Alt      Az   Elong       Time     Alt      Az   Elong";
	pheader(doc,obj,obs,"Twilight Altitude",descrip,line1,line2);
	nextDate(obsmax, dspan, obs.day);	// 'abuse' nextdate to calculate end time
	jdmax = jd(obsmax);
	while (true) {
		var jday=jd0(obscopy.year,obscopy.month,obscopy.day) + obscopy.tz/1440.0;
		if (jday >= jdmax) break;
		var rset=sunrise(obscopy,sunalt);
		// do line for current date
		doc.write(datestring(obscopy));
		for (var t=0; t<2; t++) {
			doc.write("      " + (rset[2] ? hmstring(rset[t+3],false) : "--:--"));
			if (rset[2]) {		// false if Sun never reaches specified altitude on this day
				bodies[obj].update(rset[t],obscopy);
				bodies[obj].elongupdate(rset[t],obscopy);
				doc.write("   " + fixnum(bodies[obj].alt, 5, 1));
				doc.write("  " + fixnum(bodies[obj].az, 6, 1));
				doc.write("  " + fixnum(bodies[obj].elong, 6, 1));
			}
			else doc.write("    --.-    --.-    --.-");
		}
		doc.writeln("");
		nextDate(obscopy,dstep); 
	}
	pbottom(doc,pwin,line2);
} 	// end doTwilightAltAz()


function doAltAz(obs,obj,dt) {
	// alt-az of object during one day in steps of dt minutes
	var obscopy=new Object();		// make working copy
	for (var i in obs) obscopy[i] = obs[i];
	var ostr = bodies[obj].name;
	var pwin = window.open("","pl_altitude","menubar,scrollbars,resizable");
	var doc = pwin.document;
	var line2="      Date    Time      Alt       Az     Sun    Moon";
	pheader(doc,obj,obs,"Alt-Azimuth"," ","",line2);
	var jday=jd(obscopy);
	// for each 0.5 hours do
	for (var t = 0; t < 1440/dt; t++) {
		doc.write(datestring(obscopy));
		doc.write("   " + hmstring2(obscopy.hours,obscopy.minutes,0));
		bodies[SUN].update(jday,obs);
		var sunalt = bodies[SUN].alt;
		bodies[MOON].update(jday,obs);
		var moonalt = bodies[MOON].alt;
		bodies[obj].update(jday,obs);
		doc.write("   " + fixnum(bodies[obj].alt,6,1));
		doc.write("   " + fixnum(bodies[obj].az,6,1) + "   ");
		if (sunalt > -0.833) doc.write("  Day");
		else if (sunalt > -6.0) doc.write("Civ-t");
		else if (sunalt > -12.0) doc.write("Nau-t");
		else if (sunalt > -18.0) doc.write("Ast-t");
		else doc.write(" Dark");
		if (moonalt > -0.833) doc.writeln("      Up");
		else doc.writeln("    Down");
		jday += dt/1440.0;
		adjustTime(obscopy,dt);
	}
	pbottom(doc,pwin,line2);
} 	// end doAltAz()


function doObjectData(obs,obj,dspan,dstep) {
	// list object physical data (diameter, illumination, distance, brightness)
	// eq. diameters at dist=1 AU/1km
	var obscopy=new Object(); var obsmax=new Object();
	for (var i in obs) {
		obscopy[i] = obs[i]; obsmax[i] = obs[i];
	}
	var pwin=window.open("","illumdiam","menubar,scrollbars,resizable");
	var doc=pwin.document;
	var line1="    Date      Time  Dist(Sun) Dist(Earth)    Diam.  Illum.  Brightn.";
	var line2="                         [AU]       [km]       [\"]     [%]    [mag]";
	if (obj!=MOON) 
		line2="                         [AU]       [AU]       [\"]     [%]    [mag]";
	pheader(doc,obj,obs,"Object Data"," ",line1,line2);
	nextDate(obsmax, dspan, obs.day);	// 'abuse' nextdate to calculate end time
	jdmax = jd(obsmax);
	while (true) {
		var jday=jd(obscopy);
		if (jday >= jdmax) break;
		// do line for current date
		doc.write(datestring(obscopy));
		doc.write("   " + hmstring2(obscopy.hours, obscopy.minutes,0));
		bodies[obj].update(jday,obscopy);
		var dist = bodies[obj].dist;
		var r= bodies[obj].r;
		doc.write(fixnum(r,11,3) + "    " + (obj==MOON ? fixnum(dist,7,0) : fixnum(dist,7,3)));
		if (obj<COMET) {
			doc.write(fixnum(ndiam[obj]/dist,10,1));
			doc.write(fixnum(bodies[obj].illum*100,8,1));
		} else {
			doc.write("    ----.-   ---.-");
		}
		doc.writeln(fixnum(bodies[obj].mag,9,1));
		nextDate(obscopy, dstep, obs.day); 
	}
	pbottom(doc,pwin,line2);
} 	// end doObjectData()


function doPositions(obs,obj,dspan,dstep) {
	// RA, declination, longitude, latitude and elongation
	var obscopy=new Object(); var obsmax=new Object();
	for (var i in obs) {
		obscopy[i] = obs[i]; obsmax[i] = obs[i];
	}
	var pwin=window.open("","declination","menubar,scrollbars,resizable");
	var doc=pwin.document;
	var descrip="Geocentric positions!";
	var line1="   Date       Time      RA        Dec      Alt.     Az.   Longitude  Latitude   Elongation";
	var line2="                      [h m s]    [&deg; \']      [&deg;]     [&deg;]     [&deg; \']      [&deg; \']        [&deg;]";
	pheader(doc,obj,obs,"Position",descrip,line1,line2);
	nextDate(obsmax, dspan, obs.day);	// 'abuse' nextdate to calculate end time
	jdmax = jd(obsmax);
	while (true) {
		var jday=jd(obscopy);
		if (jday >= jdmax) break;
		// do line for current date
		doc.write(datestring(obscopy));
		doc.write("   " + hmstring2(obscopy.hours, obscopy.minutes,0));
		bodies[obj].update(jday,obs);
		bodies[obj].elongupdate(jday,obs);
		var pa = bodies[obj].pa;
		doc.write("   " + hmsstring(bodies[obj].ra/15));
		doc.write("   " + anglestring(bodies[obj].dec,false,true));
		doc.write(fixnum(bodies[obj].alt,7,1));
		doc.write(fixnum(bodies[obj].az,8,1));
		doc.write("    " + anglestring(bodies[obj].eclon,true,true));
		doc.write("   " + anglestring(bodies[obj].eclat,false,true));
		doc.write(fixnum(bodies[obj].elong,8,1) + "&deg;" + (pa >= 180 ? " W" : " E"));
		doc.writeln(""); 
		nextDate(obscopy, dstep, obs.day); 
	}
	pbottom(doc,pwin,line2);
}		// end doPositions()


function doAllPlanets(obs) {
	// print data for all planets for ONE day. Positions are for local noon.
	var obscopy=new Object();
	for (var i in obs) {
		obscopy[i] = obs[i];
	}
	obscopy.hours = 0;	// set to local midnight
	obscopy.minutes = 0;
	var ndiam = [6.72, 16.68, 1, 9.36, 196.88, 165.46, 70.04, 67.0, 1, 1919.3, 716900000.0];	

	var pwin=window.open("","declination","menubar,scrollbars,resizable");
	var doc=pwin.document;
	var descrip="Date: " + datestring(obscopy);
	var line1="Object     Rise  Transit   Set       RA        Dec       Elong     Dist     Diam    Illum  Brightn";
	var line2="                                   [h m s]    [deg]      [deg]   [AU/km]     [\"]     [%]    [mag]";
	pheader(doc,-1,obs,"Data for all planets",descrip,line1,line2);
	var objects=[9,10,0,1,3,4,5,6,7];
	for (var n in objects) {
		var obj=objects[n];
		var jday=jd(obscopy);
		doc.write(bodies[obj].name);
		// rise, transit, set
		var objevents = findEvents(obj,jday,obs);
		var rise=-1; var set=-1; var objtr = -1;
		for (var i=objevents.length-2; i>0; i--) {		// scan array for relevant events
			var t=objevents[i][0]; var e = objevents[i][1];
			if (e == -1) rise = t;
			if (e == 1) set = t;
			if (e == 0) objtr = t;
		}
		doc.write("   " + dfrac2tstr(rise));
		doc.write("   " + dfrac2tstr(objtr));
		doc.write("   " + dfrac2tstr(set));
		// positions
		bodies[obj].update(jday+0.5,obs);
		bodies[obj].elongupdate(jday+0.5,obs);
		var pa = bodies[obj].pa;
		doc.write("   " + hmsstring(bodies[obj].ra/15));
		doc.write("   " + anglestring(bodies[obj].dec,false,true));
		doc.write(" " + fixnum(bodies[obj].elong,6,1) + "&deg;" + (pa >= 180 ? " W" : " E"));
		var dist = bodies[obj].dist;
		doc.write("  " + (obj==MOON ? fixnum(dist,7,0) : fixnum(dist,7,3)));
		doc.write("  " + fixnum(ndiam[obj]/dist,7,1));
		doc.write("  " + fixnum(bodies[obj].illum*100,6,1));
		doc.writeln("   " + fixnum(bodies[obj].mag,5,1));
	}
	pbottom(doc,pwin,line2);
} 	// doAllPlanets()



function doSeparation(obs,obj1,obj2,dt) {
	var obscopy=new Object();		// make working copy
	for (var i in obs) obscopy[i] = obs[i];
	var pwin = window.open("","pl_altitude","menubar,scrollbars,resizable");
	var doc = pwin.document;
	var descrip = "of " + bodies[obj1].name + " and " + bodies[obj2].name;
	var line1 = "                        Object 1       Object 2";
	var line2 = "      Date    Time     Alt     Az     Separ.     PA    Occult.   Sun";
	pheader(doc,-1,obs,"Separation",descrip,line1,line2);
	var jday=jd(obscopy);
	// for each dt minutes do
	for (var t = 0; t < 1440/dt; t++) {
		doc.write(datestring(obscopy));
		doc.write("   " + hmstring2(obscopy.hours,obscopy.minutes,0));
		bodies[SUN].update(jday,obs);
		var sunalt = bodies[SUN].alt;
		bodies[obj1].update(jday,obs);
		bodies[obj2].update(jday,obs);
		alt1=bodies[obj1].alt; alt2=bodies[obj2].alt;
		az1=bodies[obj1].az; az2=bodies[obj2].az;
		doc.write("  " + fixnum(alt1,6,1) + "  " + fixnum(az1,6,1));
		var sep = separation(az1,az2,alt1,alt2);
		doc.write("  " + fixnum(sep[0],7,2) + "&deg;" + "  " + fixnum(sep[1],6,1) + "   ");
		touchd = (ndiam[obj1]/bodies[obj1].dist + ndiam[obj2]/bodies[obj2].dist)/(2*3600);
		if (touchd > sep[0])		// check if occultation (eclipse, transit)
			doc.write("  Yes   ");
		else
			doc.write("        "); 
		if (sunalt > -0.833) doc.write("  Day");
		else if (sunalt > -6.0) doc.write("Civ-t");
		else if (sunalt > -12.0) doc.write("Nau-t");
		else if (sunalt > -18.0) doc.write("Ast-t");
		else doc.write(" Dark");
		doc.writeln("");
		jday += dt/1440.0;
		adjustTime(obscopy,dt);
	}
	pbottom(doc,pwin,line2);
} 	// end doSeparation()


function doVisible(obs,deepsky,minalt,rasort) {
// show visible stars or DSO
	var obscopy=new Object();
	for (var i in obs) {obscopy[i] = obs[i];	}
	var pwin=window.open("","visible","menubar,scrollbars,resizable");
	var doc=pwin.document;
	var descrip="Date: " + datestring(obscopy) + "  Time: " + hmstring2(obscopy.hours,obscopy.minutes,0);
	var line1="";
	var line2="Ident     Proper name         RA        Dec       Alt      Az      Rise  Transit   Set ";
	if (deepsky) 
		line2="Ident    Name          Const     RA        Dec       Alt      Az      Rise  Transit   Set ";
	pheader(doc,-1,obs,"Visible " + (deepsky?"Deep sky objects":"Stars"),descrip,line1,line2);

	var jday = jd(obscopy);
	var ord = new Array();	// for storing records of siderial time and index to catalogue
	var sid = local_sidereal(obs);
	for (var i=0; i<(deepsky?dso.length:stars.length); i++) {
		var h = parsecol(deepsky?dso[i].ra:stars[i].ra) - sid;	// negative hour angle
		if (h<0) h+=24; if (h>12) h-=24;
		ord[i] = new Array(h,i);	// 
	}
	if (rasort) isort(ord);		// sort according to transit time starting from north
	for (var k=0; k<(deepsky?dso.length:stars.length); k++) {
		i = ord[k][1];	
		var ra = (deepsky?dso[i].ra:stars[i].ra);
		var de = (deepsky?dso[i].de:stars[i].de);
		bodies[20].ra = parsecol(ra)*15;	// use User object as temporary object
		bodies[20].dec = parsecol(de);
		bodies[20].update(jday,obscopy);
		if (bodies[20].alt < minalt) continue;	// object not visible, skip it	
		if (deepsky) {
			doc.write(fixstr(dso[i].numb,9) + fixstr(dso[i].name,14) + " " + dso[i].cons);
		} 
		else {
			doc.write(fixstr(stars[i].star,6) + fixstr(stars[i].cons,4) + fixstr(stars[i].name,14))
		}
		doc.write("   " + hmsstring(bodies[20].ra/15));
		doc.write("   " + anglestring(bodies[20].dec,false,true));
		doc.write(" " + fixnum(bodies[20].alt,6,1));
		doc.write("   " + fixnum(bodies[20].az,6,1));

		var objevents = findEvents(20,jday-0.5,obs);	// center "day" on observation time
		var rise=-1; var set=-1; var objtr = -1;
		for (var j=objevents.length-2; j>0; j--) {		// scan array for relevant events
			var t=objevents[j][0]+jday-obs.tz/1440; t = t-Math.floor(t); 
			var e = objevents[j][1];
			if (e == -1) rise = t;
			if (e == 1) set = t;
			if (e == 0) objtr = t;
		}
		doc.write("    " + dfrac2tstr(rise));
		doc.write("   " + dfrac2tstr(objtr));
		doc.write("   " + dfrac2tstr(set));

		doc.writeln("");
	}

	pbottom(doc,pwin,line2);
}		// end doVisible()


function doDailyEvents(obs,obj,dspan,dstep, sunevents) {
	// List daily events (rise, set etc) for one object, optionally include solar events
	// obs is a reference variable, make a copy
	var obscopy=new Object(); var obsmax=new Object();
	for (var i in obs) {
		obscopy[i] = obs[i]; obsmax[i] = obs[i];
	}
	obscopy.hours = 0;	// set to local midnight
	obscopy.minutes = 0;
	if (dstep > 0 && dstep < 1.0) dstep = 1.0;	// time step must be at least one day
	var ostr = bodies[obj].name;
	var pwin=window.open("","daily","menubar,scrollbars,resizable");
	var doc=pwin.document;
	var line1 = "               "; var line2 = "      Date     ";
	if (obj != SUN) {
		line1 += "  |         Object        ";
		line2 += "  |   Rise  Transit   Set ";
	}
	if (obj == SUN || sunevents) {
		line1 += "  |     Morning Twilight               Sun                Evening Twilight  ";
		line2 += "  |   Astr    Naut   Civil     Rise  Transit   Set     Civil    Naut    Astr";
	}
	pheader(doc,obj,obs,"Daily Events","",line1,line2);
	nextDate(obsmax, dspan, obs.day);	// 'abuse' nextdate to calculate end time
	jdmax = jd(obsmax);
	while (true) {
		var jday=jd0(obscopy.year,obscopy.month,obscopy.day) + obscopy.tz/1440.0;
		if (jday >= jdmax) break;
		// do line for current date
		if (obj != SUN) {
			var objevents = findEvents(obj,jday,obs);
			var rise=-1; var set=-1; var objtr = -1;
			for (var i=objevents.length-2; i>0; i--) {		// scan array for relevant events
				var t=objevents[i][0]; var e = objevents[i][1];
				if (e == -1) rise = t;
				if (e == 1) set = t;
				if (e == 0) objtr = t;
			}
		}
		if (sunevents || obj == SUN) {
			var events = findEvents(SUN,jday,obs);
			srise=-1; sset=-1; ctw_b=-1; ctw_e=-1; ntw_b=-1; ntw_e=-1;  atw_b=-1; atw_e=-1; 
			for (var i=events.length-2; i>0; i--) {
				var t = events[i][0]; var e = events[i][1];
				if (e == 0) suntr = t;
				else if (e == -1) srise=t;
				else if (e == 1) sset=t;
				else if (e == -2) ctw_b = t;
				else if (e == 2) ctw_e = t;
				else if (e == -3) ntw_b = t;
				else if (e == 3) ntw_e = t;
				else if (e == -4) atw_b = t;
				else if (e == 4) atw_e = t;
			}
		}
		var dw=Math.floor(jday+1.5)-7*Math.floor((jday+1.5)/7);
		doc.write(datestring(obscopy) + "  " + dow[dw]);
		if (obj != SUN) {
			doc.write("     " + dfrac2tstr(rise) + "   " + dfrac2tstr(objtr) + "   " + dfrac2tstr(set));
		}
		if (obj == SUN || sunevents) {
			doc.write("     " + dfrac2tstr(atw_b) + "   " + dfrac2tstr(ntw_b) + "   " + dfrac2tstr(ctw_b));
			doc.write("    " + dfrac2tstr(srise) + "   " + dfrac2tstr(suntr) + "   " + dfrac2tstr(sset));
			doc.write("    " + dfrac2tstr(ctw_e) + "   " + dfrac2tstr(ntw_e) + "   " + dfrac2tstr(atw_e));
		}
		doc.writeln("");
		nextDate(obscopy,dstep); 
	}
	pbottom(doc,pwin,line2);
}		// end doDailyEvents()


;// Functions for the planets

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

;// The place, observatory definitions and daylight savings functions

// Copyright Ole Nielsen 2002-2003, Peter Hayes 1999-2001

function place(name,latitude,ns,longitude,we,zone,dss,dse) {
  this.name      = name;
  this.latitude  = latitude;
  this.ns        = ns;
  this.longitude = longitude;
  this.we        = we;
  this.zone      = zone;
  this.dss       = dss;
  this.dse       = dse;
}

// A selection of places
// Please leave Greenwich in the first entry as the default
// The second entry is my home town, I suggest you change it to yours
// is you keep a copy for your personal use.
// This database is based on Peter Hayes' original database with several places added

var atlas = new Array(
  new place("UK:Greenwich","51:28:38",0,"00:00:00",0,0,"3:5:0","10:5:0"),
  new place("NL:Rijswijk","52:02:00",0,"4:19:00",1,-60,"3:5:0","10:5:0"),
  new place("AT:Vienna","48:13:00",0,"16:22:00",1,-60,"3:5:0","10:5:0"),
  new place("AU:Melbourne","37:48:00",1,"144:58:00",1,-600,"10:5:0","03:5:0"),
  new place("AU:Perth","31:58:00",1,"115:49:00",1,-480,"10:5:0","03:5:0"),
  new place("BE:Brussels","50:50:00",0,"4:21:00",1,-60,"3:5:0","10:5:0"),
  new place("BR:Rio de Janeiro","22:54:00",1,"43:16:00",0,180,"",""),
  new place("CA:Calgary","51:03:00",0,"114:05:00",0,420,"04:1:0","10:5:0"),
  new place("CA:Halifax","44:35:00",0,"63:39:00",0,240,"04:1:0","10:5:0"),
  new place("CA:Toronto","43:39:00",0,"79:23:00",0,300,"04:1:0","10:5:0"),
  new place("CH:Zurich","47:22:40",0,"08:33:04",1,-60,"3:5:0","10:5:0"),
  new place("CL:Santiago","33:30:00",1,"70:40:00",0,240,"10:5:0","03:5:0"),
  new place("DE:Berlin","52:32:00",0,"13:25:00",1,-60,"3:5:0","10:5:0"),
  new place("DE:Frankfurt/Main","50:06:00",0,"8:41:00",1,-60,"3:5:0","10:5:0"),
  new place("DE:Hamburg","53:33:00",0,"10:00:00",1,-60,"3:5:0","10:5:0"),
  new place("DE:Munich","48:08:00",0,"11:35:00",1,-60,"3:5:0","10:5:0"),
  new place("DK:Copenhagen","55:43:00",0,"12:34:00",1,-60,"3:5:0","10:5:0"),
  new place("DK:Kolding","55:31:00",0,"9:29:00",1,-60,"3:5:0","10:5:0"),
  new place("DK:Aalborg","57:03:00",0,"9:51:00",1,-60,"3:5:0","10:5:0"),
  new place("DK:�rhus","56:10:00",0,"10:13:00",1,-60,"3:5:0","10:5:0"),
  new place("EG:Cairo","30:03:00",0,"31:15:00",1,-120,"",""),
  new place("ES:Madrid","40:25:00",0,"03:42:00",0,-60,"3:5:0","10:5:0"),
  new place("ES:Malaga","36:43:00",0,"04:25:00",0,-60,"3:5:0","10:5:0"),
  new place("ES:Las Palmas","28:08:00",0,"15:27:00",0,60,"3:5:0","10:5:0"),
  new place("FI:Helsinki","60:08:00",0,"25:00:00",1,-120,"3:5:0","10:5:0"),
  new place("FR:Bordeaux","44:50:00",0,"0:34:00",0,-60,"3:5:0","10:5:0"),
  new place("FR:Brest","48:24:00",0,"4:30:00",0,-60,"3:5:0","10:5:0"),
  new place("FR:Lille","50:38:00",0,"03:04:00",1,-60,"3:5:0","10:5:0"),
  new place("FR:Lyon","45:46:00",0,"04:50:00",1,-60,"3:5:0","10:5:0"),
  new place("FR:Marseille","43:18:00",0,"5:22:00",1,-60,"3:5:0","10:5:0"),
  new place("FR:Paris","48:48:00",0,"02:14:00",1,-60,"3:5:0","10:5:0"),
  new place("FR:Puimichel","43:58:00",0,"06:01:00",1,-60,"3:5:0","10:5:0"),
  new place("FR:Strasbourg","48:35:00",0,"7:45:00",1,-60,"3:5:0","10:5:0"),
  new place("GL:Nuuk","64:15:00",0,"51:34:00",0,180,"3:5:0","10:5:0"),
  new place("GR:Athens","38:00:00",0,"23:44:00",1,-120,"3:5:0","10:5:0"),
  new place("HK:Hong Kong","22:15:00",0,"114:11:00",1,-480,"",""),
  new place("HR:Zagreb","45:48:00",0,"15:58:00",1,-60,"3:5:0","10:5:0"),
  new place("IE:Dublin","53:19:48",0,"06:15:00",0,0,"3:5:0","10:5:0"),
  new place("IN:New Delhi","28:22:00",0,"77:13:00",1,-330,"",""),
  new place("IQ:Baghdad","33:20:00",0,"44:26:00",1,-180,"",""),
  new place("IR:Teheran","35:44:00",0,"51:30:00",1,-210,"",""),
  new place("IS:Reykjavik","64:09:00",0,"21:58:00",0,60,"3:5:0","10:5:0"),
  new place("IT:Milan","45:28:00",0,"9:12:00",1,-60,"3:5:0","10:5:0"),
  new place("IT:Palermo","38:08:00",0,"13:23:00",1,-60,"3:5:0","10:5:0"),
  new place("IT:Rome","41:53:00",0,"12:30:00",1,-60,"3:5:0","10:5:0"),
  new place("JP:Tokyo","35:70:00",0,"139:46:00",1,-540,"3:5:0","10:5:0"),
  new place("LU:Luxembourg","49:36:00",0,"6:09:00",1,-60,"3:5:0","10:5:0"),
  new place("NL:Amsterdam","52:22:23",0,"4:53:33",1,-60,"3:5:0","10:5:0"),
  new place("NL:Apeldoorn","52:13:00",0,"5:57:00",1,-60,"3:5:0","10:5:0"),
  new place("NL:Maastricht","50:51:00",0,"5:04:00",1,-60,"3:5:0","10:5:0"),
  new place("NL:Groningen","53:13:00",0,"6:33:00",1,-60,"3:5:0","10:5:0"),
  new place("NL:The Hague","52:05:00",0,"4:29:00",1,-60,"3:5:0","10:5:0"),
  new place("NL:Utrecht","52:05:10",0,"05:07:45",1,-60,"3:5:0","10:5:0"),
  new place("NO:Bergen","60:21:00",0,"5:20:00",1,-60,"3:5:0","10:5:0"),
  new place("NO:Oslo","59:56:00",0,"10:45:00",1,-60,"3:5:0","10:5:0"),
  new place("NO:Troms�","69:70:00",0,"19:00:00",1,-60,"3:5:0","10:5:0"),
  new place("NZ:Wellington","41:17:00",1,"174:47:00",1,-720,"10:5:0","03:5:0"),
  new place("PL:Warszawa","52:15:00",0,"21:00:00",1,-60,"3:5:0","10:5:0"),
  new place("PT:Faro","37:01:00",0,"7:56:00",0,0,"3:5:0","10:5:0"),
  new place("PT:Lisbon","38:44:00",0,"9:08:00",0,0,"3:5:0","10:5:0"),
  new place("PR:San Juan","18:28:00",0,"66:08:00",0,240,"04:1:0","10:5:0"),
  new place("RO:Bucharest","44:25:00",0,"26:07:00",1,-120,"3:5:0","10:5:0"),
  new place("RU:Irkutsk","52:18:00",0,"104:15:00",1,-480,"3:5:0","10:5:0"),
  new place("RU:Moscow","55:45:00",0,"37:35:00",1,-180,"3:5:0","10:5:0"),
  new place("RU:Omsk","55:00:00",0,"73:22:00",1,-360,"3:5:0","10:5:0"),
  new place("SE:Gothenburg","57:43:00",0,"11:58:00",1,-60,"3:5:0","10:5:0"),
  new place("SE:Stockholm","59:35:00",0,"18:06:00",1,-60,"3:5:0","10:5:0"),
  new place("SE:Lule�","65:35:00",0,"22:09:00",1,-60,"3:5:0","10:5:0"),
  new place("SG:Singapore","01:20:00",0,"103:50:00",1,-450,"",""),
  new place("VC:Kingstown","13:15:00",0,"61:12:00",0,240,"",""),
  new place("UK:Birmingham","52:30:00",0,"01:49:48",0,0,"3:5:0","10:5:0"),
  new place("UK:Belfast","54:34:48",0,"05:55:12",0,0,"3:5:0","10:5:0"),
  new place("UK:Cambridge","52:10:00",0,"00:06:00",0,0,"3:5:0","10:5:0"),
  new place("UK:Cardiff","51:30:00",0,"03:12:00",0,0,"3:5:0","10:5:0"),
  new place("UK:Edinburgh","55:55:48",0,"03:13:48",0,0,"3:5:0","10:5:0"),
  new place("UK:London","51:30:00",0,"00:10:12",0,0,"3:5:0","10:5:0"),
  new place("US:Anchorage","61:10:00",0,"149:53:00",0,560,"04:1:0","10:5:0"),
  new place("US:Dallas","32:48:00",0,"96:48:00",0,360,"04:1:0","10:5:0"),
  new place("US:Denver","39:45:00",0,"104:59:00",0,420,"04:1:0","10:5:0"),
  new place("US:Honolulu","21:19:00",0,"157:86:00",0,600,"04:1:0","10:5:0"),
  new place("US:Los Angeles","34:03:15",0,"118:14:28",0,480,"04:1:0","10:5:0"),
  new place("US:Miami","25:47:00",0,"80:20:00",0,300,"04:1:0","10:5:0"),
  new place("US:Minneapolis","44:58:01",0,"93:15:00",0,360,"04:1:0","10:5:0"),
  new place("US:Seattle","47:36:00",0,"122:19:00",0,480,"04:1:0","10:5:0"),
  new place("US:Washington DC","38:53:51",0,"77:00:33",0,300,"04:1:0","10:5:0"),
  new place("VC:St Vincent","13:15:00",0,"61:12:00",0,240,"",""),
  new place("ZA:Cape Town","33:56:00",1,"18:28:00",1,-120,"",""),
  new place("ZM:Lusaka","15:26:00",1,"28:20:00",1,-120,"","")
);


function observatory(place,year,month,day,hr,min,sec) {
// The observatory object holds local date and time,
// timezone correction in minutes with daylight saving if applicable,
// latitude and longitude (west is positive)
	this.name = place.name;
	this.year = year;
	this.month = month;
	this.day = day;
	this.hours = hr;
	this.minutes = min;
	this.seconds = sec;
	this.tz = place.tz;
	this.dst = false;	// is it DST?
	this.latitude = place.latitude;
	this.longitude = place.longitude;
}

// The default observatory (Greenwich noon Jan 1 2000) 
// changed by user setting place and time from menu

var observer  = new observatory(atlas[0],2000,1,1,12,0,0);

// Site name returns name and latitude / longitude as a string
function sitename() {
  var sname=observer.name;
  var latd=Math.abs(observer.latitude)+0.00001;
  var latdi=Math.floor(latd);
  sname+=((latdi < 10) ? " 0" : " ") + latdi;
  latm=60*(latd-latdi); latmi=Math.floor(latm);
  sname+=((latmi < 10) ? ":0" : ":") + latmi;
//  lats=60*(latm-latmi); latsi=Math.floor(lats);
//  sname+=((latsi < 10) ? ":0" : ":") + latsi;
  sname+=((observer.latitude >= 0) ? " N, " : " S, ");
  var longd=Math.abs(observer.longitude)+0.00001;
  var longdi=Math.floor(longd);
  sname+=((longdi < 10) ? "0" : "") + longdi;
  longm=60*(longd-longdi); longmi=Math.floor(longm);
  sname+=((longmi < 10) ? ":0" : ":") + longmi;
//  longs=60*(longm-longmi); longsi=Math.floor(longs);
//  sname+=((longsi < 10) ? ":0" : ":") + longsi;
  sname+=((observer.longitude >= 0) ? " W" : " E");
  return sname;
}	// sitename()


function checkdst(obs) {
	// Check DST is an attempt to check daylight saving, its not perfect.
	// Returns 0 or -60 that is amount to remove to get to zone time.
	// this function is now only called when selecting a place from the dropdown list. No dst check when updating the time!
	// We only know daylight saving if in the atlas
	if ((tbl.Place.selectedIndex < 0) || (tbl.Place.selectedIndex >= atlas.length))
		return 0;
	var dss=atlas[tbl.Place.selectedIndex].dss;
	var dse=atlas[tbl.Place.selectedIndex].dse;
	var ns=atlas[tbl.Place.selectedIndex].ns;
	if (dss.length==0) return 0;
	if (dse.length==0) return 0;
	// parse the daylight saving start & end dates
	var col1=dss.indexOf(":");
	var col2=dss.lastIndexOf(":");
	var col3=dss.length;
	var dssm=parseInt(dss.substring(0,col1),10);
	var dssw=parseInt(dss.substring(col1+1,col2),10);
	var dssd=parseInt(dss.substring(col2+1,col3),10);
	col1=dse.indexOf(":");
	col2=dse.lastIndexOf(":");
	col3=dse.length;
	var dsem=parseInt(dse.substring(0,col1),10);
	var dsew=parseInt(dse.substring(col1+1,col2),10);
	var dsed=parseInt(dse.substring(col2+1,col3),10);
	// Length of months
	// year,month,day and day of week
	var jdt=jd0(obs.year,obs.month,obs.day);
	var ymd=jdtocd(jdt);
	// first day of month - we need to know day of week
	var fymd=jdtocd(jdt-ymd[2]+1);
	// look for daylight saving / summertime changes
	// first the simple month checks
	// Test for the northern hemisphere
	if (ns==0) {
		if ((ymd[1]>dssm) && (ymd[1]<dsem)) return -60;
		if ((ymd[1]<dssm) || (ymd[1]>dsem)) return 0;
	} 
	else{
		// Southern hemisphere, New years day is summer.
		if ((ymd[1]>dssm) || (ymd[1]<dsem)) return -60;
		if ((ymd[1]<dssm) && (ymd[1]>dsem)) return 0;
	}
	// check if we are in month of change over
	if (ymd[1]==dssm) { // month of start of summer time
		// date of change over
		var ddd=dssd-fymd[3]+1;
		ddd=ddd+7*dssw;
		while (ddd>month_length[ymd[1]-1]) ddd-=7;
		if (ymd[2]<ddd) return 0;
		// assume its past the change time, its impossible
		// to know if the change has occured.
		return -60;
	} 
	if (ymd[1]==dsem) { // month of end of summer time
		// date of change over
		var ddd=dsed-fymd[3]+1;
		ddd=ddd+7*dsew;
		while (ddd>month_length[ymd[1]-1]) ddd-=7;
		if (ymd[2]<ddd) return -60;
		// see comment above for start time
		return 0;
	}
	return 0;
}	// checkdst()


function jd(obs) {
// The Julian date at observer time
  var tz = obs.tz || 0;

	var j = jd0(obs.year,obs.month,obs.day);
	j+=(obs.hours+((obs.minutes+tz)/60.0)+(obs.seconds/3600.0))/24;
	return j;
}	// jd()


function local_sidereal(obs) {
// sidereal time in hours for observer
	var res=g_sidereal(obs.year,obs.month,obs.day);
	res+=1.00273790935*(obs.hours+(obs.minutes+obs.tz+(obs.seconds/60.0))/60.0);
	res-=obs.longitude/15.0;
	while (res < 0) res+=24.0;
	while (res > 24) res-=24.0;
	return res;
}

;// Various date and time functions

// Copyright Peter Hayes 1999-2001, Ole Nielsen 2002-2004


// must be updated using leapyear() if year changed
var month_length=new Array(31,28,31,30,31,30,31,31,30,31,30,31);
var dow = ["Sun","Mon","Tue","Wed","Thu","Fri","Sat"];


function leapyear(year) {
  var leap=false;
  if (year % 4 == 0) leap = true;
  if (year % 100 == 0 ) leap = false;
  if (year % 400 == 0) leap = true;
  return leap;
}


function jd0(year,month,day) {
// The Julian date at 0 hours(*) UT at Greenwich
// (*) or actual UT time if day comprises time as fraction
  var y  = year;
  var m = month;
  if (m < 3) {m += 12; y -= 1};
  var a = Math.floor(y/100);
  var b = 2-a+Math.floor(a/4);
  var j = Math.floor(365.25*(y+4716))+Math.floor(30.6001*(m+1))+day+b-1524.5;
  return j;
}	// jd0()


function jdtocd(jd) {
// The calendar date from julian date, see Meeus p. 63
// Returns year, month, day, day of week, hours, minutes, seconds
  var Z=Math.floor(jd+0.5);
  var F=jd+0.5-Z;
  if (Z < 2299161) {
    var A=Z;
  } else {
    var alpha=Math.floor((Z-1867216.25)/36524.25);
    var A=Z+1+alpha-Math.floor(alpha/4);
  }
  var B=A+1524;
  var C=Math.floor((B-122.1)/365.25);
  var D=Math.floor(365.25*C);
  var E=Math.floor((B-D)/30.6001);
  var d=B-D-Math.floor(30.6001*E)+F;
  if (E < 14) {
    var month=E-1;
  } else {
    var month=E-13;
  }
  if ( month>2) {
    var year=C-4716;
  } else {
    var year=C-4715;
  }
  var day=Math.floor(d);
  var h=(d-day)*24;
  var hours=Math.floor(h);
  var m=(h-hours)*60;
  var minutes=Math.floor(m);
  var seconds=Math.round((m-minutes)*60);
  if (seconds >= 60) {
    minutes=minutes+1;
    seconds=seconds-60;
  }
  if (minutes >= 60) {
    hours=hours+1;
    minutes=0;
  }
  var dw=Math.floor(jd+1.5)-7*Math.floor((jd+1.5)/7);
  return new Array(year,month,day,dw,hours,minutes,seconds);  
}	// jdtocd()


function g_sidereal(year,month,day) {
// sidereal time in hours for Greenwich
  var T=(jd0(year,month,day)-2451545.0)/36525;
  var res=100.46061837+T*(36000.770053608+T*(0.000387933-T/38710000.0));
  return rev(res)/15.0;
}

;// Extensions to the Math routines - Trig routines in degrees

// Copyright Peter Hayes 1999-2001, Ole Nielsen 2003-2004

var DEG2RAD = Math.PI/180.0;
var RAD2DEG = 180.0/Math.PI;

function rev(angle) 	{return angle-Math.floor(angle/360.0)*360.0;}		// 0<=a<360
function rev2(angle)	{var a = rev(angle); return (a>=180 ? a-360.0 : a);}	// -180<=a<180
function sind(angle) 	{return Math.sin(angle*DEG2RAD);}
function cosd(angle) 	{return Math.cos(angle*DEG2RAD);}
function tand(angle) 	{return Math.tan(angle*DEG2RAD);}
function asind(c) 		{return RAD2DEG*Math.asin(c);}
function acosd(c) 		{return RAD2DEG*Math.acos(c);}
function atand(c) 		{return RAD2DEG*Math.atan(c);}
function atan2d(y,x) 	{return RAD2DEG*Math.atan2(y,x);}

function log10(x) 		{return Math.LOG10E*Math.log(x);}

function sqr(x)			{return x*x;}
function cbrt(x)		{return Math.pow(x,1/3.0);}

function SGN(x) 		{ return (x<0)?-1:+1; }

;// SUN and MOON
// Alternative version of Sun position based on Schlyter's method

// Copyright Ole Nielsen 2002-2004
// Please read copyright notice in astrotools2.html source

// 'Meeus' means "Astronomical Algorithms", 2nd ed. by Jean Meeus

// ecliptic position of the Sun relative to Earth (basically simplified version of planetxyz calc)
function sunxyz(jday) {
	// return x,y,z ecliptic coordinates, distance, true longitude
	// days counted from 1999 Dec 31.0 UT
	var d=jday-2451543.5;
	var w = 282.9404 + 4.70935E-5 * d;		// argument of perihelion
	var e = 0.016709 - 1.151E-9 * d; 
	var M = rev(356.0470 + 0.9856002585 * d); // mean anomaly
	var E = M + e*RAD2DEG * sind(M) * ( 1.0 + e * cosd(M) );
	var xv = cosd(E) - e;
	var yv = Math.sqrt(1.0 - e*e) * sind(E);
	var v = atan2d( yv, xv );		// true anomaly
	var r = Math.sqrt( xv*xv + yv*yv );	// distance
	var lonsun = rev(v + w);	// true longitude
	var xs = r * cosd(lonsun);		// rectangular coordinates, zs = 0 for sun 
	var ys = r * sind(lonsun);
	return new Array(xs,ys,0,r,lonsun,0);
}

function SunAlt(jday,obs) {
// return alt, az, time angle, ra, dec, ecl. long. and lat=0, illum=1, 0, dist, brightness 
	var sdat=sunxyz(jday);
	var ecl = 23.4393 - 3.563E-7 * (jday-2451543.5); 
	var xe = sdat[0];
	var ye = sdat[1]*cosd(ecl);
	var ze = sdat[1]*sind(ecl);
	var ra = rev(atan2d(ye,xe));
	var dec = atan2d(ze, Math.sqrt(xe*xe+ye*ye));
	var topo=radec2aa(ra,dec,jday,obs);
	return new Array( topo[0] , topo[1] , topo[2] , ra , dec , sdat[4], 0, 1, 0, sdat[3], -26.74 );
}


// Sun rise and set times (if twilight==-0.833) or desired twilight time. Return julian days
function sunrise(obs,twilight) {
	// obs is a reference variable make a copy
	var obscopy=new Object();
	for (var i in obs) obscopy[i] = obs[i];
	obscopy.hours=12;
	obscopy.minutes=0;
	obscopy.seconds=0;
	var riseset=new Array(0.0, 0.0, false, 0.0, 0.0);
	var lst=local_sidereal(obscopy);
	var jday=jd(obscopy);
	var radec = SunAlt(jday, obscopy);
	var ra = radec[3]; var dec = radec[4];
	var UTsun=12.0+ra/15.0-lst;
	if (UTsun < 0.0) UTsun+=24.0;
	if (UTsun > 24.0) UTsun-=24.0;
	var cosLHA = (sind(twilight)-sind(obs.latitude)*sind(dec)) / (cosd(obs.latitude)*cosd(dec));
	// Check for midnight sun and midday night. "riseset[2]" false if no rise and set found
	riseset[2]=false;
	if (cosLHA <= 1.0 && cosLHA >= -1.0) {
		// rise/set times allowing for not today.
		riseset[2]=true;
		var lha=acosd(cosLHA)/15.0;
		if ((UTsun-lha) < 0.0) {
			var rtime=24.0+UTsun-lha;
		} else {
			var rtime=UTsun-lha;
		}
		riseset[0]=jday+rtime/24.0-0.5;
		if ((UTsun+lha) > 24.0) {
			var stime=UTsun+lha-24.0;
		} else {
			var stime=UTsun+lha;
			riseset[4]=stime;
		}
		riseset[1]=jday+stime/24.0-0.5;
		// riseset[3] and [4] are times in UT hours
		riseset[3]=rtime;
		riseset[4]=stime;
	}
	return(riseset);
}


function MoonPos(jday,obs) {
// MoonPos calculates the Moon position and distance, based on Meeus chapter 47
// and the illuminated percentage from Meeus equations 48.4 and 48.1
// OPN: This version of MoonPos calculates the position to a precision of about 2' or so
// All T^2, T^3 and T^4 terms skipped. NB: Time is not taken from obs but from jday (julian day)
// Returns alt, az, hour angle, ra, dec (geocentr!), eclip. long and lat (geocentr!), 
// illumination, distance, brightness and phase angle 
	var T=(jday-2451545.0)/36525;
	// Moons mean longitude L'
	var LP=rev(218.3164477+481267.88123421*T);
	// Moons mean elongation
	var D=rev(297.8501921+445267.1114034*T);
	// Suns mean anomaly
	var M=rev(357.5291092+35999.0502909*T);
	// Moons mean anomaly M'
	var MP=rev(134.9633964+477198.8675055*T);
	// Moons argument of latitude
	var F=rev(93.2720950+483202.0175233*T);
	// The "further arguments" A1, A2 and A3  and the term E have been ignored
	// Sum of most significant terms from table 45.A and 45.B (terms less than 0.004 deg / 40 km dropped)
	var Sl = 6288774*sind(MP) + 1274027*sind(2*D-MP) + 658314*sind(2*D) + 213618*sind(2*MP) -
		185116*sind(M) - 114332*sind(2*F) + 58793*sind(2*D-2*MP) + 57066*sind(2*D-M-MP) +
		53322*sind(2*D+MP) + 45758*sind(2*D-M) - 40923*sind(M-MP) - 34720*sind(D) - 
		30383*sind(M+MP) + 15327*sind(2*D-2*F) - 12528*sind(MP+2*F) + 10980*sind(MP-2*F) +
		10675*sind(4*D-MP) + 10034*sind(3*MP) + 8548*sind(4*D-2*MP) - 7888*sind(2*D+M-MP) -
		6766*sind(2*D+M) - 5163*sind(D-MP) + 4987*sind(D+M) + 4036*sind(2*D-M+MP);
	var Sb = 5128122*sind(F) + 280602*sind(MP+F) + 277602*sind(MP-F) + 173237*sind(2*D-F) +
		55413*sind(2*D-MP+F) + 46271*sind(2*D-MP-F) + 32573*sind(2*D+F) + 17198*sind(2*MP+F) +
		9266*sind(2*D+MP-F) + 8822*sind(2*MP-F) + 8216*sind(2*D-M-F) + 4324*sind(2*D-2*MP-F) +
		4200*sind(2*D+MP+F);
	var Sr = (-20905355)*cosd(MP) - 3699111*cosd(2*D-MP) - 2955968*cosd(2*D) - 569925*cosd(2*MP) +
		246158*cosd(2*D-2*MP) - 152138*cosd(2*D-M-MP) - 170733*cosd(2*D+MP) - 204586*cosd(2*D-M) -
		129620*cosd(M-MP) + 108743*cosd(D) + 104755*cosd(M+MP) + 79661*cosd(MP-2*F) + 48888*cosd(M);
	// geocentric longitude, latitude
	var mglong = rev(LP+Sl/1000000.0);
	var mglat = Sb/1000000.0;
	// Obliquity of Ecliptic
	var obl = 23.4393-3.563E-7*(jday-2451543.5);
	var ra = rev(atan2d(sind(mglong)*cosd(obl)-tand(mglat)*sind(obl),cosd(mglong)));
	var dec = asind(sind(mglat)*cosd(obl)+cosd(mglat)*sind(obl)*sind(mglong));
	var moondat=radec2aa(ra,dec,jday,obs);
	// phase angle (48.4)
	var pa = Math.abs(180.0 - D - 6.289*sind(MP) + 2.100*sind(M) - 1.274*sind(2*D-MP) - 
		0.658*sind(2*D) - 0.214*sind(2*MP) - 0.11*sind(D));
	var k = (1+cosd(pa))/2;
	var mr = Math.round(385000.56+Sr/1000.0);
	var h = moondat[0];
	// correct for parallax, equatorial horizontal parallax, see Meeus p. 337
	h -= asind(6378.14/mr)*cosd(h);
	// brightness, use Paul Schlyter's formula (based on common phase law for Moon)
	var sdat=sunxyz(jday);
	var r = sdat[3];	// Earth (= Moon) distance to Sun in AU
	var R = mr/149598000;	// Moon distance to Earth in AU
	var mag = 0.23 + 5*log10(r*R) + 0.026 * pa + 4.0E-9 * pa*pa*pa*pa
	return new Array(h,moondat[1],moondat[2],ra,dec,mglong, mglat, k, r, mr, mag);
}	// Moonpos()


;// Utility functions

// Copyright Ole Nielsen 2002-2004, Peter Hayes 1999-2001


function datestring(obs) {
// datestring provides a locale independent format
  var datestr = "";  datestr += obs.year;
  datestr += ((obs.month < 10) ? ":0" : ":") + obs.month;
  datestr += ((obs.day < 10) ? ":0" : ":") + obs.day;
  return datestr;
}		// end datestring()


function datestring2(year,month,day) {
  var datestr = "";  datestr += year;
  datestr += ((month < 10) ? ":0" : ":") + month;
  datestr += ((day < 10) ? ":0" : ":") + day;
  return datestr;
}		// end datestring2()


function adjustTime(obs,amount) {
// update date and time, amount is in minutes (may be negative)
// added 2004
    month_length[1] = leapyear(obs.year) ? 29 : 28;
	if (amount<0) {
		amount = Math.abs(amount);
		obs.minutes -= amount%60; amount = Math.floor(amount/60.0);
		obs.hours -= amount%24; amount = Math.floor(amount/24.0);
		obs.day -= amount;
		if (obs.minutes < 0) {
			obs.minutes += 60;
			obs.hours -= 1;
		}
		if (obs.hours < 0) {
			obs.hours += 24;
			obs.day -= 1;
		}
		while (obs.day < 1) {
			obs.day += month_length[obs.month-2+(obs.month==1?12:0)];
			obs.month -= 1; 
			if (obs.month == 0) {
				obs.year -= 1;
				obs.month = 12;
				month_length[1] = (leapyear(obs.year) ? 29 : 28);
			}
		}
	}
	else {
		obs.minutes+=amount%60; amount=Math.floor(amount/60.0);
		obs.hours+=amount%24; amount=Math.floor(amount/24.0);
		obs.day+=amount;
		if (obs.minutes > 59) {
			obs.minutes -= 60;
			obs.hours += 1;
		}
		if (obs.hours > 23) {
			obs.hours -= 24;
			obs.day += 1;
		}
		while (obs.day > month_length[obs.month-1]) {
			obs.day -= month_length[obs.month-1];
			obs.month += 1; 
			if (obs.month == 13) {
				obs.year += 1;
				obs.month = 1;
				month_length[1] = (leapyear(obs.year) ? 29 : 28);
			}
		}
	}
}		// end adjustTime()


function hmsstring(t) {
// the caller must add a leading + if required.
  var hours = Math.abs(t);
  var minutes = 60.0*(hours-Math.floor(hours));
  hours=Math.floor(hours);
  var seconds = Math.round(60.0*(minutes-Math.floor(minutes)));
  minutes=Math.floor(minutes);
  if (seconds >= 60) { minutes+=1; seconds-=60; }
  if (minutes >= 60) { hour+=1; minutes-=60; }
  if (hours >= 24) { hours-=24; }
  var hmsstr=(t < 0) ? "-" : "";
  hmsstr=((hours < 10) ? "0" : "" )+hours;
  hmsstr+=((minutes < 10) ? ":0" : ":" )+minutes;
  hmsstr+=((seconds < 10) ? ":0" : ":" )+seconds;
  return hmsstr;
}		// end hmsstring()


function hmstring(t,plus) {
// hmstring converts hours to a string (+/-)hours:minutes, used for relative time (TZ)
  var hours = Math.abs(t);
  var minutes = Math.round(60.0*(hours-Math.floor(hours)));
  hours=Math.floor(hours);
  if (minutes >= 60) { hours+=1; minutes-=60; }	// minutes could be 60 due to rounding
  if (hours >= 24) { hours-=24; }
  var hmstr = (t < 0) ? "-" : (plus?"+":"");
  hmstr += ((hours < 10) ? "0" : "" )+hours;
  hmstr += ((minutes < 10) ? ":0" : ":" )+minutes;
  return hmstr;
}		// end hmstring()


function hmstring2(hours,minutes,seconds) {
// hmstring2 returns time as a string HH:MM (added 2004.01.02), seconds needed for rounding
	if (seconds>=30) minutes++;
	if (minutes>=60) {hours++; minutes=0;}
	var timestr = ((hours < 10) ? "0" : "") + hours;	
	timestr += ((minutes < 10) ? ":0" : ":") + minutes;
	return timestr;
}		// end hmstring2()


function dmsstring(d) {
// dmsstring converts lat/long angle to unsigned string d:m:s
  var deg = Math.abs(d);
  var minutes = 60.0*(deg-Math.floor(deg));
  deg=Math.floor(deg);
  var seconds = Math.round(60.0*(minutes-Math.floor(minutes)));
  minutes=Math.floor(minutes);
  if (seconds >= 60) { minutes+=1; seconds-=60; }
  if (minutes >= 60) { deg+=1; minutes-=60; }
  hmsstr=((deg < 10) ? "0" : "" )+deg;
  hmsstr+=((minutes < 10) ? ":0" : ":" )+minutes;
  hmsstr+=((seconds < 10) ? ":0" : ":" )+seconds;
  return hmsstr;
}		// end dmsstring()


function dmstring(d) {
// dmstring converts lat/long angle to unsigned string d:m
  var deg = Math.abs(d);
  var minutes = 60.0*(deg-Math.floor(deg));
  deg=Math.floor(deg);
  var seconds = Math.round(60.0*(minutes-Math.floor(minutes)));
  minutes=Math.floor(minutes);
  if (seconds >= 30) { minutes+=1; }
  if (minutes >= 60) { deg+=1; minutes-=60; }
  hmstr=((deg < 10) ? "0" : "" )+deg;
	hmstr+=((minutes < 10) ? ":0" : ":" )+minutes;
  return hmstr;
} // end dmstring()


function anglestring(a,circle,arcmin) {
	// Return angle as degrees:minutes. 'circle' is true for range between 0 and 360 
	// and false for -90 to +90, if 'arcmin' use deg and arcmin symbols
  var ar=Math.round(a*60)/60;
  var deg=Math.abs(ar);
  var min=Math.round(60.0*(deg-Math.floor(deg)));
  if (min >= 60) { deg+=1; min=0; }
  var anglestr="";
  if (!circle) anglestr+=(ar < 0 ? "-" : "+");
  if (circle) anglestr+=((Math.floor(deg) < 100) ? "0" : "" );
  anglestr+=((Math.floor(deg) < 10) ? "0" : "" )+Math.floor(deg);
  if (arcmin) anglestr+=((min < 10) ? "&deg;0" : "&deg;")+(min)+"' ";
	else anglestr+=((min < 10) ? ":0" : ":" )+(min);
  return anglestr;
} // end anglestring()


function fixnum(n,l,d) {
	// convert float n to right adjusted string of length l with d digits after decimal point.
	// the sign always requires one character, allow for that in l!
	var m = 1;
	for (var i=0; i<d; i++) m*=10;
	var n1 = Math.round(Math.abs(n)*m);
	var nint = Math.floor(n1/m);
	var nfract = (n1 - m*nint) + ""; // force conversion to string
	while (nfract.length < d) nfract = "0" + nfract;
	var str = (n<0 ? "-" : " ") + nint;
	if (d > 0) str = str + "." + nfract;
	while (str.length<l) str = " " + str;
	return str;
} // end fixnum()


function fixstr(str,l) {
	// returns left-adjusted string of length l, pad with spaces or truncate as necessary
	if (str.length > l) return str.substring(0,l); 
	while (str.length < l) {
		str += " ";
	}
	return str;
}		// end fixstr()


function parsecol(str) {
	// parsecol converts deg:min:sec or hr:min:sec to a number
  var col1=str.indexOf(":");
  var col2=str.lastIndexOf(":");
  if (col1 < 0) return parseInt(str);
  if (str.substring(0,1) == "-") {
    var res=parseInt(str.substring(1,col1),10);
  } else {
    var res=parseInt(str.substring(0,col1),10);
  }
  if (col2 > col1) {
    res+=(parseInt(str.substring(col1+1,col2),10)/60.0) +
         (parseInt(str.substring(col2+1,str.length),10)/3600.0);
  } else {
    res+=(parseInt(str.substring(col1+1,str.length),10)/60.0);
  }
  if (str.substring(0,1) == "-") {
    return -res;
  } else {
    return res;
  }
}	// end parsecol()


function interpol(n,y1,y2,y3) {
	// interpolate y (Meeus 3.3)
	var a = y2-y1;
	var b = y3-y2;
	var c = b-a;
	return y2+(n/2)*(a+b+n*c);
}


function nzero(y1,y2,y3) {
	// Calculate value of interpolation factor for which y=zero. n0 should be within [-1:1] 
	// Meeus formula (3.7)	
	var a = y2-y1; var b = y3-y2; var c = b-a;
	var n0 = 0;
	do {
		dn0 = -(2*y2 + n0*(a+b+c*n0)) / (a+b+2*c*n0);
		n0 += dn0;
	} while (Math.abs(dn0) > 0.0001);
	return n0;
}		// end nzero()


function nextrem(y1,y2,y3) {
	// Calculate value of interpolation factor for which y reaches extremum (-1<n<1);
	var a = y2-y1; var b = y3-y2; var c = b-a;
	var nm = -(a+b)/(2*c);	// (3.5)
	return nm;
}		// end nextrem();


function isort(arr) {
// Sort 2D array in ascending order on first column of each element using insertion sort 
	for (var c=0; c<arr.length-1; c++) {	
		var tmp = arr[c+1]; var a = c;
		while (a>=0 && arr[a][0]>tmp[0]) {
			arr[a+1] = arr[a];
			a--;
		}
		arr[a+1] = tmp;
	}
}	// end isort()


