# spec/integration/plaac_web_spec.rb
require 'spec_helper'
require 'nokogiri'
require 'rack/test'

RSpec.describe 'PLAAC web workflow', type: :request do

  def run_web_workflow(fasta_name)
    fasta_path = File.expand_path("../../../cli/example/#{fasta_name}", __dir__)

    # Step 1: POST file to /find with required parameters
    file = Rack::Test::UploadedFile.new(fasta_path, 'text/plain')

    post '/find', file: file, len: 60, alpha: 100

    # puts "=== POST /find response ==="
    # puts "Status: #{last_response.status}"
    # puts "Headers: #{last_response.headers.inspect}"
    # puts "Body (first 500 chars): #{last_response.body[0..500]}"

    # Step 2: Follow redirects to results page (if any)
    while last_response.redirect?
      location = last_response.headers['Location']
      # puts "Redirected to: #{location}"
      get location
    end

    # puts "=== Final response after redirects ==="
    # puts "Status: #{last_response.status}"
    # puts "Headers: #{last_response.headers.inspect}"
    # puts "Body (first 500 chars): #{last_response.body[0..500]}"

    # Step 3: Parse HTML for TSV link
    html = Nokogiri::HTML(last_response.body)
    tsv_link_element = html.at_css('a[href*="/candidates/"][href$="/tsv"]')
    expect(tsv_link_element).not_to be_nil

    tsv_link = tsv_link_element['href']
    puts "TSV link found: #{tsv_link}"

    # Step 4: GET the TSV URL (simulating user clicking the link)
    get tsv_link
    expect(last_response).to be_ok
    expect(last_response.content_type).to include('text/tab-separated-values')

    yield html
  end

  def check_visualization_images(html)
    images = html.css('img[src]')
    expect(images.length).to eq(5)
    image_urls = images.map { |img| img['src'] }
    expect(image_urls).to all(match(%r{\A/visualize/[^/]+/(?:strippng|images/[1-4]\.png)\z}))

    expect(image_urls).to include(
                            a_string_matching(%r{/strippng\z}),
                            a_string_matching(%r{/images/1\.png\z}),
                            a_string_matching(%r{/images/2\.png\z}),
                            a_string_matching(%r{/images/3\.png\z}),
                            a_string_matching(%r{/images/4\.png\z})
                          )

    image_urls.each do |image_url|
      get image_url
      expect(last_response).to be_ok
      expected_type =
        if image_url.end_with?('/strippng')
          'application/png'
        else
          'image/png'
        end
      expect(last_response.content_type).to eq(expected_type)
      puts "Found #{image_url}, #{last_response.content_type}"
    end
  end

  def visualize_all_candidates(html)
    token = html.at_css('input[name="token"]')['value']
    
    selected = html.css('#candidate_list input[name="selected[]"]')
                 .map { |input| input['value'] }
    
    # puts "Visualize token: #{token}"
    # puts "Selected candidates: #{selected.inspect}"

    header 'Cookie', 'plaac_finder_picklist=[1,2,3,4]'    
    post '/visualize',
       token: token,
       'selected[]' => selected
    
    # puts "=== POST /visualize response ==="
    # puts "Status: #{last_response.status}"
    # puts "Headers: #{last_response.headers.inspect}"
    # puts "Body (first 1000 chars): #{last_response.body[0..1000]}"
    
    while last_response.redirect?
      location = last_response.headers['Location']
      # puts "Redirected to: #{location}"
      get location
      # puts "GET status: #{last_response.status}"
      # puts "GET body (first 1000 chars): #{last_response.body[0..1000]}"
    end
    
    expect(last_response).to be_ok
    Nokogiri::HTML(last_response.body)
  end  
  
  it 'uploads a file, runs analysis, and downloads the correct TSV via app URLs' do
    run_web_workflow('MOT3.fasta') do
      # Optional debug: print first few lines of TSV
      # puts "TSV content (first 500 chars):"
      # puts last_response.body[0..500]

      # Step 5: Compare downloaded TSV content to gold file
      # Note this strips out version number, since that will change
      gold_path = File.expand_path('../../../cli/example/MOT3-candidates-plaac-1.1.0.tsv', __dir__)

      gold_content = File.read(gold_path)

      actual = last_response.body.gsub(
        /^## plaac_version=.*$/,
        '## plaac_version=<version>;'
      )

      expected = gold_content.gsub(
        /^## plaac_version=.*$/,
        '## plaac_version=<version>;'
      )
      expect(actual).to eq(expected)
    end
  end

  it 'generates visualization PNGs and checks them for four_classic_prions.fasta' do
    run_web_workflow('four_classic_prions.fasta') do |html|
      visualize_html = visualize_all_candidates(html)
      check_visualization_images(visualize_html)
    end
  end

  it 'generates visualization PNGs and checks them for four_classic_prions_tabs.fasta' do
    run_web_workflow('four_classic_prions_tabs.fasta') do |html|
      visualize_html = visualize_all_candidates(html)
      check_visualization_images(visualize_html)
    end
  end
end
