path="/Users/simon/Desktop/WIKIPEDIA_PROCESS"
`ls #{path}/HMM_SOL`.split("\n").select { |i| i.include?("dat") }.each { |i|
  file=File.new("#{path}/HMM_SOL/#{i}", 'r')
  file_out=File.new("#{path}/HMM_SOL_new/#{i}", 'w')
  file_out.write(file.readline)
  str=file.readline
  n=str.split("=")[1].to_i
  file_out.write(str)
  file_out.write(file.readline)
  n.times { 
    file_out.write(file.readline)
  }
  n.times {
    file_out.write(file.readline)    
  }
  file_out.write("CR\n")
  file_out.write(file.read)
  file.close
  file_out.close
}

class String
  def renorm(scale=2)
    n=self.length
    str_new=""
    (n/scale).times { |i|
      str_new << self[i*scale]
    }
    str_new
  end
end

`rm -rf #{path}/RENORM/*`

filename="God_12states"
1.upto(10) { |scale|
  print "***** Doing #{scale} renormalization\n"
  test_file="#{path}/HMM_SOL_new/#{filename}_hmm.dat"
  str=`./hmm -g #{50000*scale} #{test_file}`.split("\n")[1].renorm(scale)
  file=File.new("#{path}/RENORM/#{filename}_#{scale}renorm", 'w'); file.write("#{str.length}\n#{str}"); file.close
  print `./parallel_glue.rb 8 10 #{path}/RENORM/#{filename}_#{scale}renorm`
  while(`ls #{path}/RENORM/#{filename}_#{scale}renorm_OUT_*states 2>/dev/null`.split("\n").length < 1) do
    sleep(1)
  end
  print "Found renormalized system for #{scale}\n"
  print "Going on...\n"
}

