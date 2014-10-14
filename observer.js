// The place, observatory definitions and daylight savings functions

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
  new place("DK:Århus","56:10:00",0,"10:13:00",1,-60,"3:5:0","10:5:0"),
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
  new place("NO:Tromsø","69:70:00",0,"19:00:00",1,-60,"3:5:0","10:5:0"),
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
  new place("SE:Luleå","65:35:00",0,"22:09:00",1,-60,"3:5:0","10:5:0"),
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

