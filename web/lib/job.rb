class Job

  STDOUT_FILE="stdout.log"
  STDERR_FILE="stderr.log"

  WORKING_DIRECTORY_BASE = "/tmp/clusty"


  attr_accessor :command, :owner, :status, :started_at, :id, :exit_status, :ended_at, :script_name
  attr_reader :token
  # -----------------------------------

  def initialize(token=nil)
    @token = token || UUID.new.generate
  end

  def stdout
    @stdout ||= File.read(File.join(working_directory,STDOUT_FILE))
    @stdout
  end

  def stderr
    @stderr ||= File.read(File.join(working_directory,STDERR_FILE))
    @stderr
  end

  def params=(hash)
    create_working_directory
    File.open(File.join(working_directory,'params.yml'), 'w') do |file|
      file.write YAML.dump(hash)
    end
    nil
  end

  def params
    @params ||= YAML.load(IO.read(File.join(working_directory,'params.yml')))
    @params
  rescue Errno::ENOENT
    @params
  end

  def update_params(key,val)
    h = params
    h[key] = val
    self.params = h
    @params = nil
    nil
  end

  def working_directory
    File.join(WORKING_DIRECTORY_BASE, token)
  end
  def create_working_directory(dir=working_directory)
    unless File.directory?(dir)
      system "mkdir -p #{dir}"
    end
  end
  def add_file(filename,dir=working_directory)
    create_working_directory(dir)
    `cp #{filename} #{dir}`
  end

  def save()
    true
  end

  def submit()
    Job.perform(self)
  end

  # -----------------------------------
  class << self
    @@log = Logger.new(File.expand_path('../../logs/plaac.log',__FILE__))

    private

    def log
      @@log
    end

    def log_block(level,tag,msg)
      unless msg.nil?
        msg.each_line do |line|
          log.send(level, "%s %s" % [tag, line.strip])
        end
      end
    end
  end

  def self.enqueue(job)
    Job.perform(job)
  end

  def self.first(opts={})
    Job.new(opts[:token])
  end

  def self.perform(job)
    job.status = "running"
    job.started_at = DateTime.now
    save_ok = job.save
    log.debug("running.save #{job.inspect}")


    FileUtils.mkpath(job.working_directory)
    script_file = File.join(job.working_directory, job.script_name || "run.sh")
    File.open(script_file, 'w') do |file|
      file.puts job.command
    end
    FileUtils.chmod(0755, script_file)
    cmd = %Q[ cd #{job.working_directory} && sh #{script_file} 1>stdout.log 2>stderr.log ]
    log.info("LOCAL RUN %s %s [%s] %s" % [job.id, job.token, job.owner, job.command.inspect])
    %x[ #{cmd} ]
    log.info("RESULT-EXIT %s pid=%d exit=%d" % [job.id, $?.pid, $?.exitstatus])

    job.exit_status = $?.exitstatus.to_i
    job.status = ($?.exitstatus == 0) ? "success" : "failure"
    job.ended_at = DateTime.now

    log_block(:info, "RESULT-STDOUT", job.stdout)
    stderr = job.stderr
    log_block(:warn, "RESULT-STDERR", stderr) unless stderr.nil?

    save_ok = job.save
    log.debug("finished.save #{job.inspect}")

    job
  rescue Exception => e
    puts e
  end
end
