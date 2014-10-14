module.exports = function(grunt) {

  grunt.initConfig({
    pkg: grunt.file.readJSON('package.json'),
    concat: {
      options: {
        separator: ';',
        stripBanners: true,
      },
      lib: {
        src : [
          'makelist.js',
          'planets.js',
          'observer.js',
          'datetime.js',
          'math.js',
          'sunmoon.js',
          'util.js'
        ],
        dest: 'blastro.js'
      }
    }
  });

  grunt.loadNpmTasks('grunt-contrib-concat');


  ////////////////////////////////////////////////////
  // Main grunt tasks
  ////////////////////////////////////////////////////

  grunt.registerTask('build', ['concat']);

};
