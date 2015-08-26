class ProteinsController < ApplicationController

  def index
    @protein = Protein.all
  end

  def new
    @protein = Protein.new
  end

  def show
  end

  def create






    require 'net/http'
    require 'uri'

      rna_codons = {"UUU" => 'F', "UUC" => 'F', "UUA" => 'L', "UUG" => 'L',
                    "UCU" => 'S', "UCC" => 'S', "UCA" => 'S', "UCG" => 'S',
                    "UAU" => 'Y', "UAC" => 'Y', "UAA" => '-', "UAG" => '-',
                    "UGU" => 'C', "UGC" => 'C', "UGA" => '-', "UGG" => 'W',
                    "CUU" => 'L', "CUC" => 'L', "CUA" => 'L', "CUG" => 'L',
                    "CCU" => 'P', "CCC" => 'P', "CCA" => 'P', "CCG" => 'P',
                    "CAU" => 'H', "CAC" => 'H', "CAA" => 'Q', "CAG" => 'Q',
                    "CGU" => 'R', "CGC" => 'R', "CGA" => 'R', "CGG" => 'R',
                    "AUU" => 'I', "AUC" => 'I', "AUA" => 'I', "AUG" => 'M',
                    "ACU" => 'T', "ACC" => 'T', "ACA" => 'T', "ACG" => 'T',
                    "AAU" => 'N', "AAC" => 'N', "AAA" => 'K', "AAG" => 'K',
                    "AGU" => 'S', "AGC" => 'S', "AGA" => 'R', "AGG" => 'R',
                    "GUU" => 'V', "GUC" => 'V', "GUA" => 'V', "GUG" => 'V',
                    "GCU" => 'A', "GCC" => 'A', "GCA" => 'A', "GCG" => 'A',
                    "GAU" => 'D', "GAC" => 'D', "GAA" => 'E', "GAG" => 'E',
                    "GGU" => 'G', "GGC" => 'G', "GGA" => 'G', "GGG" => 'G'}

      # Fetches SEQUENCE from GenBank
      # Sets ID (temporary. Will have user input the ID or name in the simple_form on proteins#new)
      accession = 296334

      # need to get the sequence
      # then convert from DNA to RNA,
      #


      # # Defines the URI for eSearch
      # uri = URI('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi')

      # # Defines params for dynamic GET request for eSearch
      # params = { :db => "nucleotide", :id => accession, :sort => "relevance" }

      # uri.query = URI.encode_www_form(params)

      # response = Net::HTTP.get_response(uri)
      # seq = response.body if response.is_a?(Net::HTTPSuccess)

        # seq4 = seq3.join.delete("\n")

      # Defines the URI
      # Defines params for dynamic GET request
      # bobby recommends httparty
      def accession(accession)
        uri = URI('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi')
        params = { :db => "protein", :id => accession, :rettype => "fasta", :retmode => "text"}
        uri.query = URI.encode_www_form(params)
        response = Net::HTTP.get_response(uri)
        seq = response.body if response.is_a?(Net::HTTPSuccess)
      end

      def convertDNA(seq)
        seqArray = []
        seq.split(/\n/).each do |i|
          seqArray << /\A((?!>.*).)*\z/.match(i)
        end
        seq = seqArray.join
      end

      def convertyThyUr
        # Replaces Ts with Us
        seq.gsub!("T","U")
      end

      # Calls the appropriate residue from the hash
      def convertAA(seq)
        i = 0
        residuesSeq = []
        while i < (seq.length / 3)
          residuesSeq << rna_codons[seq[i..(i + 2)]]
          i += 3
        end
      end

      def mw(residuesSeq)
        weight = []
        residuesSeq.each do |i|
          weight << mw[i].to_f
        end
      end


      puts "The molecular weight of the protein is #{weight.inject(:+)}"

      # puts seq


      puts residuesSeq.join


  end

  def update
  end
end
