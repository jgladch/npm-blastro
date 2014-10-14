// Extensions to the Math routines - Trig routines in degrees

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

