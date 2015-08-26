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

dna_complement = {'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C'}

# MW is for mass upon addition (water is already gone)
mw = {'A' => 71.03711, 'C' => 103.00919, 'D' => 115.02694, 'E' => 129.04259,
      'F' => 147.06841, 'G' => 57.02146, 'H' => 137.05891, 'I' => 113.08406,
      'K' => 128.09496, 'L' => 113.08406, 'M' => 131.04049, 'N' => 114.04293,
      'P' => 97.05276, 'Q' => 128.05858, 'R' => 156.10111, 'S' => 87.03203,
      'T' => 101.04768, 'V' => 99.06841, 'W' => 186.07931, 'Y' => 163.06333,
      # denotes carboxyamidomethylation, ~ denotes phosphorylation, @ denotes acetylation
      '#' => 57.02146, '~' => 79.96633, '@' => 42.01056}

# pKa dictionaries
acid_chains      = {'C' => 8.14, 'D' => 3.71, 'E' => 4.15}
base_chains      = {'H' => 6.04, 'K' => 10.67, 'R' => 12.10}
side_chains      = {'C' => 8.14, 'D' => 3.71, 'E' => 4.15, 'H' => 6.04, 'K' => 10.67, 'R' => 12.10}
amino_groups     = {'A' => 9.71, 'C' => 10.28, 'D' => 9.66, 'E' => 9.58,
                    'F' => 9.09, 'G' => 9.58, 'H' => 9.09, 'I' => 9.60,
                    'K' => 9.16, 'L' => 9.58, 'M' => 9.08, 'N' => 8.76,
                    'P' => 10.47, 'Q' => 9.00, 'R' => 9.00, 'S' => 9.05,
                    'T' => 8.96, 'V' => 9.52, 'W' => 9.34, 'Y' => 9.04}
carboxyl_groups  = {'A' => 2.33, 'C' => 1.91, 'D' => 1.95, 'E' => 2.16,
                    'F' => 2.18, 'G' => 2.34, 'H' => 1.70, 'I' => 2.26,
                    'K' => 2.15, 'L' => 2.32, 'M' => 2.16, 'N' => 2.16,
                    'P' => 1.95, 'Q' => 2.18, 'R' => 2.03, 'S' => 2.13,
                    'T' => 2.20, 'V' => 2.27, 'W' => 2.38, 'Y' => 2.24}

# Fetches SEQUENCE from GenBank
# Sets ID (temporary. Will have user input the ID or name in the simple_form on proteins#new)
accession = 296334


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

# Calculates the molecular weight of the peptide



# puts weight.inject { |sum, n| sum + n }

# weight = []
# i = 0
# while i < (residuesSeq.length)
#   puts (mw[i])
#   i += 1
# end





# puts residuesSeq.join















# need to match any line that doesn't start with \>

# splits by newlines and joins
# puts (seq.split(/\n/)).join

# ack = seq.split(/\n/)
# puts ack.class
# puts ack

# right now I have everything as a string. I need to join everything that doesn't start with \A>

# split everything that matches the regex. the regex matches all except lines that start with '>'. then join.

# /^(>.*)$/
# ^([^>].*)$ possibly this

# seq = seq.split(/\n/)
# puts /^([^>].*)$/.match(seq)

  # seq2 = seq.split(/\n/)
  # seq3 = []
  # seq2.each do |i|
  #   seq3 << /^((?!>.*).)*$/.match(i)
  # end


# seq4 = []
# seq4 = /^((?!>.*).)*$/.match(seq.split(/\n/))
# puts seq4

# seq3 = []
# seq3 = seq.scan(/^((?!>.*).)*$/)
# puts seq3.join

# seq2 = seq.split(/\n/)
# puts /^((?!>.*).)*$/.match(seq2)



# puts /\A>.\z/.match(seq.split(/\n/))


# puts /(G.*\n)/.match(seq)

# puts seq

# puts seq
# puts seq.class

# puts /\n/.match(seq)
# puts /^([^<].*)$/.match(seq.split("\n"))

# puts /^[A-Z]+$/.match(seq)


# puts /(^[^>].*$)/.match(seq)

# puts   /\A[^\>]/.match(seq)

# if
#   # \A may also be ^ for start of line
#   /\A[^\>]/.match(seq)



# puts seq
# puts seq.length
# puts seq.length / 3
# puts seq[0..2]
# puts seq[2..4]


















# def open(url)
#   Net::HTTP.get(URI.parse(url))
# end

# page_content = open('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=296334&rettype=fasta&retmode=text')
# puts page_content





# sequence = []
# sequence << "#{seq.scan(/.{3}/)}"

# puts sequence.class

# sequence.each do |i|
#   puts rna_codons["#{i}"]
# end

# sequence.each do |i|
#   puts rna_codons["#{i}"]
# end




# sequence = []

# seq.each_char do |i|
#   if i == "T"
#     i = "U"
#     sequence << "#{i}"
#   else
#     sequence << i
#   end
#   # sequence.join
# end

# puts sequence.join

# def residues(seq)
  # puts seq.scan(/.{3}/)
# end
