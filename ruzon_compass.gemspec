Gem::Specification.new do |s|
  s.name        = 'ruzon_compass'
  s.version     = '0.0.0'
  s.date        = '2012-01-21'
  s.summary     = "Ruzon Compass Edge Detector"
  s.description = "A simple hello world gem"
  s.authors     = ["Jordan Phillips"]
  s.email       = 'jordan@bonanza.com'
  s.files		  = Dir.glob('lib/**/*.rb') +
            	    Dir.glob('ext/**/*.{c,h,rb}')
  s.extensions  = ['ext/ruzon_compass/extconf.rb']

  s.homepage    =
    'http://rubygems.org/gems/ruzon_compass'
end