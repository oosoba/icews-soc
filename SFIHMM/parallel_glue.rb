#!/usr/bin/ruby

n_cores=ARGV[0].to_i
n_iterations=ARGV[1].to_i
filename=ARGV[2]

if ARGV[3] != nil then
  start=ARGV[3].to_i
  finish=ARGV[4].to_i
else
  start=1
  finish=20
end

# n_cores=4
# n_iterations=10
# filename="test_data"

`rm *MULTI* 2>/dev/null`
n_cores.times { |i|
  `rm #{filename}_MULTI#{i}* 2>/dev/null`
  `cp #{filename} #{filename}_MULTI#{i}`
}

current_best=-1e12
list=[]
start_time=Time.now
start.upto(finish) { |states|
  print "Doing #{states} states.\n"
  start_time_state=Time.now
  n_cores.times { |i|
    `./hmm -f #{n_iterations} #{states} #{filename}_MULTI#{i} > /dev/null &`
  }
  begin
    while(`ls #{filename}_MULTI*_OUT_#{states}states 2>/dev/null`.split("\n").length < n_cores) do
      sleep(1)
    end
    # print "ALL HMMs for state #{states} finished running -- #{`ls #{filename}_MULTI*_OUT_#{states}states 2>/dev/null`.split("\n").length}\n"
    # print `ls #{filename}_MULTI*_OUT_#{states}states 2>/dev/null`
  rescue Interrupt
    print "Killing background processes\n"
    `killall hmm`
    exit
  end
  print "Done running\n"
  best_logl=Array.new(n_cores) { |i|
    file=File.new("#{filename}_MULTI#{i}_OUT_#{states}states")
    best=file.readline.split("=")[1].to_f
    file.close
    [i,best,best-states*states,states]
  }.sort { |i,j|
    i[1] <=> j[1]
  }[-1]
  n_cores.times { |i|
    if i != best_logl[0] then
      `rm #{filename}_MULTI#{i}_OUT_#{states}states`
    end
  }
  list << best_logl
  state_best_aic=best_logl[1]-states*states
  print "Current AIC: #{state_best_aic}\n"
  print "Took #{(Time.now-start_time_state)/(60)} minutes.\n"
  if state_best_aic > current_best then
    current_best=state_best_aic
  else
    break
  end
}
print "Total time: #{(Time.now-start_time)/60} minutes\n"
found=list.collect { |i| 
  i[3]
}
save_states=found[-2]
save_core=list.select { |i| i[3] == save_states }[0][0]
save_logl=list.select { |i| i[3] == save_states }[0][1]

found.each { |i|
  if i != save_states then
    `rm -rf #{filename}_MULTI*_OUT_#{i}states`
  end
}
`mv #{filename}_MULTI#{save_core}_OUT_#{save_states}states #{filename}_OUT_#{save_states}states`
`rm #{filename}_MULTI*`

print "Found! Best fit is #{save_states} states; log-l=#{save_logl}; filename #{filename}_OUT_#{save_states}states\n"
if (save_states == finish) then
  print "Unusual: the maximum number of states to try is #{finish}; you've exceeded default state size; you may wish to re-run with the command"
  print "./parallel_glue.rb #{n_cores} #{n_iterations} #{filename} #{finish*5}\n"
end
