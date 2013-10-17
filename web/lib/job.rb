class Job
  STDOUT_FILE="stdout.log"
  STDERR_FILE="stderr.log"

  WORKING_DIRECTORY_BASE = "/tmp/clusty"
  attr_accessor :command, :owner, :status, :started_at, :id, :exit_status, :ended_at
  attr_reader :token
  # -----------------------------------

  def initialize(token=nil)
    @token = token || UUID.new.generate
  end

  def stdout
    @stdout ||= File.read(File.join(@working_directory,STDOUT_FILE)) rescue nil
    @stdout
  end

  def stderr
    @stderr ||= File.read(File.join(@working_directory,STDERR_FILE)) rescue nil
    @stderr
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
    private
    def log(action,msg)
      @@logger ||= Logger.new('/tmp/clusty.log')
      $stderr.puts "%s %06s %s" % [DateTime.now.strftime("%Y-%m-%d %H:%M:%S,%L"), action, msg]
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
    log("DEBUG","running.save #{job.inspect}")


    FileUtils.mkpath(job.working_directory)
    script_file = File.join(job.working_directory,"run.sh")
    File.open(script_file, 'w') do |file|
      file.puts job.command
    end
    FileUtils.chmod(0755, script_file)
    cmd = %Q[ cd #{job.working_directory} && ./run.sh 1>stdout.log 2>stderr.log ]
    log("LOCAL RUN","%s %s [%s] %s" % [job.id, job.token, job.owner, job.command.inspect])
    %x[ #{cmd} ]
    log("RESULT","%s pid=%d exit=%d" % [job.id, $?.pid, $?.exitstatus])

    job.exit_status = $?.exitstatus.to_i
    job.status = ($?.exitstatus == 0) ? "success" : "failure"
    job.ended_at = DateTime.now
    save_ok = job.save
    log("DEBUG"," finished.save #{job.inspect}")

    job
  rescue Exception => e
    puts e
  end
end
