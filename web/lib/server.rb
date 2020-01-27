# Encoding: utf-8
require 'bundler/setup'
require 'sinatra/base'
require 'sinatra/cookies'
require 'fileutils'
require 'logger'
require 'bio'
require 'yaml'
require 'uuid'
require 'better_errors'
require 'benchmark'

$:.push File.expand_path("..",__FILE__)
require 'assets'
require 'job'

if RUBY_VERSION =~ /1.9/ # assuming you're running Ruby ~1.9
  Encoding.default_external = Encoding::UTF_8
  Encoding.default_internal = Encoding::UTF_8
end

# ensure required directories
%w[tmp logs].each do |dir|
  dir = File.expand_path("../../"+dir,__FILE__)
  FileUtils.mkdir(dir) unless File.directory?(dir)
end

ENVIRONMENT = ENV["RAILS_ENV"] || "development"
$version = `git log | head -n 1`.split(/\s/)[1][0,10] rescue "?"

class Server < Sinatra::Base
  helpers Sinatra::Cookies
  register Assets

  enable :sessions
  set :session_secret, 'marshmallow'
  set :bind, '0.0.0.0'

  configure :development do
    use BetterErrors::Middleware
    BetterErrors.application_root = File.dirname(__FILE__)
  end

  $config = YAML.load_file(File.expand_path("../../config.yaml",__FILE__))

  @@base = File.expand_path('../..',__FILE__)

  @@log = Logger.new(@@base+"/logs/plaac.log")

  @@input_fasta="input_sequence.fasta"
  @@plaac_candidates="plaac_candidates.tsv"
  @@plaac_candidates_details="plaac_candidates_details.tsv"
  @@plaac_candidates_details_pdf="plaac_candidates.pdf"

  @@background_meta = "bg_freqs/background_metadata.txt"

  bg_freq_filename = ->(id){ id.nil? ? nil : "bg_freqs_#{id}.txt" }
  bg_freq_local_filename = ->(id){ id.nil? ? nil : "bg_freqs/#{bg_freq_filename[id]}" }

  before do
    # Strip the last / from the path
    request.env['PATH_INFO'].gsub!(/\/$/, '')
  end

  #--------------------------------------------------------------------
  # Pages

  get '/' do
    @bg_freqs = read_bg_freqs()
    haml :index
  end

  get '/details' do
    haml :details
  end

  #--------------------------------------------------------------------
  # Actions

  # upload a fasta file for processing
  post '/find' do
    @@log.info("POST /run")
    @bg_freqs = read_bg_freqs()
    # validation
    unless params[:len] && params[:len].to_i > 0
      @error = "Core length must be a positive integer"
      @@log.info(@error)
      return haml(:index)
    end

    core_len = params[:len].to_i
    alpha = params[:alpha].to_f
    alpha = 100 if alpha > 100
    alpha = alpha / 100.0 rescue 0.0
    alpha 

    job = Job.new
    job.script_name = "find_candidates.sh"
    working_directory = job.working_directory
    output_filename = File.join(working_directory,@@plaac_candidates)
    filename = File.join(working_directory,@@input_fasta)

    job.add_file(path("bin/plaac_plot.r"), working_directory) # r library
    job.add_file(path("bin/plaac_plot_util.r"), working_directory) # r library
    job.add_file(path("bin/plaac.jar"), working_directory) # binary

    # Create input file on disk.
    if params[:file] &&
        (tmpfile = params[:file][:tempfile]) &&
        (name = params[:file][:filename])
      # Accept File (uploaded)
      write_file(filename, tmpfile)
    elsif params[:sequence] && !params[:sequence].empty?
      # Accept String (pasted)
      #@@log.debug "raw sequence \"#{params[:sequence].inspect}\""
      write_string(filename, add_gene_id(params[:sequence]))
    else
      @error = "No sequence."
      @@log.info(@error)
      return haml(:index)
    end

    @@log.debug "validating fasta #{filename}"
    valid_fasta, msg = valid_fasta_file(filename)
    @@log.debug "valid?: #{valid_fasta}, #{msg}"
    if !valid_fasta
      @error = msg
      return haml(:index)
    end

    gene_count = `egrep '^>' #{filename} |wc -l`.to_i

    @@log.info "gene_count: #{gene_count}"

    bg_freqs = read_bg_freqs() # ignore warning, bg_freqs used by filename lambdas
    bgfreq_filename = bg_freq_filename[params[:bg]]
    bg = params[:bg] if @bg_freqs.keys.include? params[:bg]

    if !bgfreq_filename.nil? && params[:bg] != "input_sequence"
      job.add_file(path(bg_freq_local_filename[params[:bg]]), working_directory)
      bgfreq_opt = " -B #{bgfreq_filename} "
    end
    job.params = {
      core_len: core_len,
      alpha: alpha,
      bgfreq_filename: bgfreq_filename,
      bg: bg,
      bgfreq_opt: bgfreq_opt,
    }


    job.command = <<-COMMAND.gsub(/\n/,'').squeeze(" ")
      java -jar #{working_directory}/plaac.jar -i #{filename} -c #{core_len} -a #{alpha} #{bgfreq_opt} > #{output_filename};
      touch #{working_directory}/results_ready
    COMMAND

    begin
      job.submit()
    rescue Exception => e
      @@log.error(e)
    end

    # If we only have one gene, skip to visualization
    if gene_count == 1
      params['token'] = job.token
      #cookies[:plaac_finder_picklist] = "1"
      response.set_cookie(:plaac_finder_picklist,
                          value: '1',
                          domain: nil,
                          path: '/')
      job.update_params(:picklist_ids, [1])

      display_plaac_candidates()
    else
      redirect "#{$config[:apppath]}/candidates/#{job.token}"
    end
  end

  get '/candidates/:token/tsv' do
    @job = Job.first(:token => params[:token])
    send_file File.join(@job.working_directory,"plaac_candidates.tsv"),
      :type => 'text/tab-separated-values',
      :filename => "plaac_candidates.tsv"
  end

  get '/candidates/:token' do
    @@log.info "candidates #{params[:token]}"
    @job = Job.first(:token => params[:token])
    # TODO: check that @job existed
    results_ready_file = File.join(@job.working_directory, "results_ready")
    plaac_candidates_file = File.join(@job.working_directory, @@plaac_candidates)


    # get header, and produce a hash of column symbol to array index.
    if File.exist?(results_ready_file)
      @@log.info "results ready"
      lines = File.open(plaac_candidates_file).each_line
      lines = lines.reject{|line| line =~ /^##/}
      header_line = lines.first
      if header_line.nil?
        sleep 2
        redirect "#{$config[:apppath]}/candidates/#{@job.token}?bounce"
        return
      end
      @header = header_line.split(/\t/)
      @idx = Hash[@header.each_with_index.map{|column, i| [column.to_sym, i]}]

      # TODO: REVIEW: We should move this sorting into the java.
      sorted_flag_file = File.join(@job.working_directory, ".candiates.sorted")

      unless File.exist?(sorted_flag_file)
        i = 10
        until File.writable?(plaac_candidates_file) && File.exist?(plaac_candidates_file) || i == 0
          i -= 1
          sleep 1
        end
        @@log.info "sorting #{plaac_candidates_file}"
        # get rid of first line first, as it is header, already grabbed above
        @output = lines.drop(1).map{|line| line.split(/\t/)}
       

        # Sort the candidates: COREscore descending, LLR descending (/w NaNs at the end).
        @output = @output.sort_by{|line| 
          corescore = line[@idx[:COREscore]]
          llr = line[@idx[:LLR]]
          [ (corescore=="NaN" ? 1 : 0), -1*corescore.to_f,
                  (llr=="NaN" ? 1 : 0), -1*llr.to_f
          ]
        }
      
        
        # save re-sorted file to make sure they are in the right order for next step
        File.open(plaac_candidates_file, 'w') do |s|
          s.puts @header.join("\t")
          @output.each do |candidate|
            s.puts candidate.join("\t")
          end
        end

        FileUtils.touch sorted_flag_file
      end

      @candidates, @picklist = load_candidates(@job, plaac_candidates_file)
    else
      redirect "#{$config[:apppath]}/candidates/#{@job.token}"
      return
    end

    @@log.info "render!"
    haml :candidates
  end

  post '/visualize' do
    display_plaac_candidates()
  end

  get '/visualize/:token/pdf' do
    @job = Job.first(:token => params[:token])
    send_file File.join(@job.working_directory,"plaac_details.pdf"),
      :type => 'application/pdf',
      :filename => "plaac_details.pdf"
  end

  get '/visualize/:token/strippdf' do
    @job = Job.first(:token => params[:token])
    send_file File.join(@job.working_directory,"plaac_strip.pdf"),
      :type => 'application/pdf',
      :filename => "plaac_strip.pdf"
  end

  get '/visualize/:token/strippng' do
    @job = Job.first(:token => params[:token])
    send_file File.join(@job.working_directory,"plaac_strip.png"),
      :type => 'application/png',
      :filename => "plaac_strip.png"
  end

  get '/visualize/:token/images/:id.png' do
    @job = Job.first(:token => params[:token])
    if @job.nil?
      @@log.warn("unknown image token #{params[:token]}")
      ""
    else
      image_filename="plaac_details_#{"%05d" % [params[:id].to_i]}.png"
      file = File.join(@job.working_directory,image_filename)
      send_file file,
        :type => 'image/png',
        :filename => image_filename
    end
  end

  get '/visualize/:token/tsv' do
    @job = Job.first(:token => params[:token])
    send_file File.join(@job.working_directory,@@plaac_candidates_details),
      :type => 'text/tab-separated-values',
      :filename => "plaac_details.tsv"
  end

  get '/visualize/:token' do
    @job = Job.first(:token => params[:token])
    # p @job.params
    candidates_filename = File.join(@job.working_directory,"plaac_candidates_selected.tsv")
    @output = File.open(candidates_filename).readlines.map{|line| line.split(/\t/)}


    # TODO: this is stupid, we should be using the plaac_finder_picklist to random access the file, not load the whole thing.
    plaac_candidates_file = File.join(@job.working_directory, @@plaac_candidates)
    lines = IO.readlines(plaac_candidates_file)
    lines = lines.reject{|line| line =~ /^##/}
    header_line = lines.first

    @header = header_line.try(:split, /\t/)
    @candidates = []
    unless @header.nil?
      @idx = Hash[@header.each_with_index.map{|column, i| [column.to_sym, i]}]
      all_candidates, @picklist = load_candidates(@job, plaac_candidates_file)

      selected_rows = @job.params[:picklist_ids]

      all_candidates.each_with_index do |candidate, i|
        @candidates << candidate[0].split(/\t/) if selected_rows.include?(i)
      end
    end

    haml :visualize
  end

  #--------------------------------------------------------------------
  # Helpers

  # an absolute path for filename relative to @@base
  def path(filename)
    File.expand_path("#{filename}",@@base)
  end

  # writes filename with the contents from string
  def write_string(filename,string)
    File.open(filename,'wb+') do |f|
      f.write scrub(string)
    end
  end

  # writes filename with the contents from iostream
  def write_file(filename,iostream)
    File.open(filename,'wb+') do |f|
      while (blk = iostream.read(65536))
          f.write scrub(blk)
      end
    end
  end

  def scrub(str)
    str.gsub!(/ /,' ') # remove non-breaking space
    str
  end

  def time(comment='',&block)
    a = Time.now.to_f
    block.call
    b = Time.now.to_f
    p [comment,b-a]
  end

  def load_candidates(job,filename)
    @dataset = []
    @picklist = []

    page = params[:page] || 1
    filter = nil
    filter = /#{params[:filter]}/ unless [nil,''].include? params[:filter].nil?
    filtered = []
    picklist_ids = []

    cookie_picklist = cookies[:plaac_finder_picklist]
    if cookie_picklist.nil? 
      if !job.params[:picklist_ids].nil?
        picklist_ids = job.params[:picklist_ids]
      else
        @@log.warn "no picklist!"
      end
    else
      begin
        cookie_picklist.gsub!(/[^0-9\,]/,'') # clean input
        picklist_ids = cookie_picklist.split(/,/).map(&:to_i)
        job.update_params(:picklist_ids, picklist_ids)
      rescue Error => e
        @@log.warn("Could not load cookie list: #{e}")
        picklist_ids = []
      end
    end
    @@log.debug "picklist_ids.size = #{picklist_ids.size}"

    if filter
      # filter with grep
      all_lines = IO.readlines(filename)
      all_lines = all_lines.reject{|line| line =~ /^##/}
      all_lines = all_lines.each_with_index.map{|line,i| [line,i]}
      @picklist = all_lines.select{|pair| picklist_ids.include? pair[1] }
      filtered = all_lines.select{|pair| pair[0] =~ filter}

      @total = filtered.size - 2
    else
      header_lines = 2
      @total = %x[ wc -l #{filename} ].strip.split(/ /).first.to_i - header_lines
    end

    @per_page = 50

    @pages = (@total.to_f / @per_page).ceil

    @page = page.to_i || 1
    @page = 1 if @page > @pages

    start = (@page - 1) * @per_page
    start = 0 if start < 0

    if filter
      @dataset = filtered[start, @per_page]
    else
      # very fast if you want everything
      file = File.open(filename)
      all_lines = file.lines.each_with_index.map{|line,i| [line,i]}

      @picklist = all_lines.select{|pair| picklist_ids.include? pair[1] }

      (start).times { all_lines.next } # skip skip skip
      @dataset = all_lines.take(@per_page)
    end

    @@log.debug "load_candidates size:[%s, %s]" % [@dataset.size, @picklist.size]
    [@dataset, @picklist]
  end

  def add_gene_id(sequence)
    raw_sequences = []
    raw_sequence = ""
    lines = sequence.strip.lines
    #@@log.debug "sequence_lines: #{lines.to_a.size}"
    lines.each do |line|
      if line.strip == ""
        unless raw_sequence.empty?
          raw_sequences << raw_sequence
          raw_sequence = ""
        end
      else
        raw_sequence += line
      end
    end
    unless raw_sequence.empty?
      raw_sequences << raw_sequence
    end

    raw_sequences = raw_sequences.map{|sequence| sequence.strip}

    sequences = []
    raw_sequences.each do |sequence|
      if sequence[0] == '>'
        sequences << sequence
      else
        #sequences << ">gene#{UUID.generate[0..7]}\r\n"+sequence
        sequences << ">sequence\r\n"+sequence
      end
    end
    final_sequence = sequences.join("\n\n")
    #@@log.debug "final_sequence: #{final_sequence.inspect}"
    final_sequence
  end

  def page_url(page)
    "?"+ { :filter => params[:filter], :page => page }.map{|e| e.join("=")}.join("&")
  end

  private

  def local(cmd)
    @@log.info "plaac running system(#{cmd})"
    system(cmd)
  end

  def valid_fasta_file(filename)
    invalid_fasta_rex = /([^A-Z*-])/
    File.open(filename).each_with_index do |line,i|
      next if line =~ /^>/
      return [false, "Found invalid amino acid on line #{i}. '#{$1}'"] if line.strip =~ invalid_fasta_rex
    end
    return [true, nil]
  end

  def display_plaac_candidates
    @@log.info "display_plaac_candidates #{params.inspect}"
    @job = Job.first(:token => params["token"])
    @@log.debug @job
    @candidates = []

    plaac_candidates_file = File.join(@job.working_directory, @@plaac_candidates)

    core_len = @job.params[:core_len]
    alpha = @job.params[:alpha]

    # passed in from previous page
    bgfreq_filename = @job.params[:bgfreq_filename]
    bg = @job.params[:bg]

    # rewrite plaac_candidates, with selections
    selected = [params[:selected]].flatten.compact.map(&:to_i)
    @job.update_params(:selected, selected)

    candidates_filename = File.join(@job.working_directory,"plaac_candidates_selected.tsv")
    File.open(candidates_filename,'w') do |candidates_file|
      candidates, picklist = load_candidates(@job, plaac_candidates_file)
      picklist.each do |columns,i|
        name = columns.split(/\t/)[0]
        # TODO: ultimately the second "name" should be replaced by the a user-specified "common" gene name
        candidates_file.puts name + "\t" + name
        @candidates << name
      end
    end
    local("chmod ug+rw #{candidates_filename}")

    input_fasta = @@input_fasta
    output_filename = @@plaac_candidates_details

    job = Job.new(@job.token)
    job.script_name = "visualize.sh"
    working_directory = job.working_directory

    # generate details
    if !bgfreq_filename.nil? && bg != "input_sequence"
      bg_freq_filename = ->(id){ id.nil? ? nil : "bg_freqs_#{id}.txt" }
      bg_freq_local_filename = ->(id){ id.nil? ? nil : "bg_freqs/#{bg_freq_filename[id]}" }
      job.add_file(path(bg_freq_local_filename[bg]), working_directory)
      bgfreq_opt = " -B #{bgfreq_filename} "
    end

    job.command = <<-COMMAND
    sh -c "cd #{@job.working_directory} &&
      java -jar ./plaac.jar -i #{input_fasta} -c #{core_len} -a #{alpha} -p #{candidates_filename} #{bgfreq_opt} > #{output_filename} &&
      ./plaac_plot.r plaac_candidates_details.tsv plaac_details.pdf &&
      ./plaac_plot.r plaac_candidates_details.tsv plaac_details.png &&
      ./plaac_plot.r plaac_candidates_details.tsv plaac_strip.pdf -c &&
      ./plaac_plot.r plaac_candidates_details.tsv plaac_strip.png -c ;
      touch details_ready
    "
    COMMAND
    #job.working_directory = @job.working_directory
    job.owner = session[:login]
    if job.save
      @@log.info "ok"
    else
      job.errors.each do |error|
        @@log.error(error)
      end
    end

    job.submit()

    redirect "#{$config[:apppath]}/visualize/#{job.token}"
  end

  def read_bg_freqs
    bg_freqs = {}
    skipfirst = true
    File.open(@@background_meta).each do |line|
      if skipfirst
        skipfirst = false
        next
      end
      (uniprot_version, id_name, species_name, fasta_url, download_date, aa_freq_filename) = line.split(/\t+/)
      bg_freqs[id_name] = species_name
    end
    bg_freqs
  end

  run! if __FILE__ == $0
end
