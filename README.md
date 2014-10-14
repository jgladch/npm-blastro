# Blastro #

> A Node module for generating planetary ephemeris data for a given date and location. This module repurposes code originally published by Ole Nielsen in 2002 and Peter Hayes in 1999.

## Summary ##
  > Blastro will generate ephemeris data for planetary bodies for any given date and location. This data can be used for astronomical and astrological pursuits.

## Contributing ##
> This is a barebones implementation of one of Ole Nielsen's ephemeris functions. I plan to continue to develop this module and aim to create the standard for javascript ephemerides.

## Usage ##
> Currently, this module exports a single function (returnAllPlanets). The single argument to this function is an object in the following format:

<pre><code>obs = {
    day: 14,                //integer
    dst: false,             //daylight savings time, boolean
    hours: 6,               //integer
    latitude: "42:43:38",   //string in this format
    longitude: "82:43:00",  //string in this format
    minutes: 50,            //integer
    month: 11,              //integer
    name: "Location",       //string
    seconds: 0,             //integer
    tz: 0,                  //timezone correction from GMT
    year: 2014              //integer
  };</code></pre>
