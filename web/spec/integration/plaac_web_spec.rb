# spec/integration/plaac_web_spec.rb
require 'spec_helper'
require 'nokogiri'
require 'rack/test'

RSpec.describe 'PLAAC web workflow', type: :request do
  it 'uploads a file, runs analysis, and downloads the correct TSV via app URLs' do
    fasta_path = File.expand_path('../../../cli/example/MOT3.fasta', __dir__)
    gold_path  = File.expand_path('../../../cli/example/MOT3-candidates-plaac-1.1.0.tsv', __dir__)

    # Step 1: POST file to /find with required parameters
    file = Rack::Test::UploadedFile.new(fasta_path, 'text/plain')

    post '/find', file: file, len: 60, alpha: 100

    puts "=== POST /find response ==="
    puts "Status: #{last_response.status}"
    puts "Headers: #{last_response.headers.inspect}"
    puts "Body (first 500 chars): #{last_response.body[0..500]}"

    # Step 2: Follow redirects to results page (if any)
    while last_response.redirect?
      location = last_response.headers['Location']
      puts "Redirected to: #{location}"
      get location
    end

    puts "=== Final response after redirects ==="
    puts "Status: #{last_response.status}"
    puts "Headers: #{last_response.headers.inspect}"
    puts "Body (first 500 chars): #{last_response.body[0..500]}"

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

    # Optional debug: print first few lines of TSV
    puts "TSV content (first 500 chars):"
    puts last_response.body[0..500]

    # Step 5: Compare downloaded TSV content to gold file
    # Note this strips out version number, since that will change
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
