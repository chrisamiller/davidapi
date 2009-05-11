#!/usr/bin/env ruby

# Author:  Chris Miller
# Contact: chrisamiller@gmail.com
#
# This program is free software; you can redistribute it
# and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more
# details.
#
# The full text for the General Public License is here:
# http://www.gnu.org/licenses/gpl.html
#


require 'rubygems'
require 'mechanize'
require 'choice'

#------------------------------------------------------------------------------
# Takes an appropriately formatted API URL and retrieves the plaintext
# table returned by DAVD
#
def retrieveTable(url)
  if url.length > 2048
    raise "URL length is limited to 2048 characters. Use a smaller gene list and try again"
  end

  #retrieve the results page from the API
  puts "calling the DAVID API with the following URL:"
  puts url
  puts ""

#  WWW::Mechanize.html_parser = Nokogiri
  mech = WWW::Mechanize.new { |agent|
    agent.user_agent_alias = 'Mac Safari'
  }

  #get the page
  page = mech.get(url)

  #parse out the script lines containing relevant variables
  script = []
  page.search("script").each do |input|
    script << input.content
  end

  script = script[0]

  rowids = ""
  annot = ""
  action = ""

  script.each{|line|

    if line =~ /rowids/
      rowids = line.split(/=/)[1].gsub(/\"/,"").gsub(/\;/,"")
      rowids.strip!
    elsif line =~/annot\.value/
      annot = line.split(/=/)[1].gsub(/\"/,"").gsub(/\;/,"")
      annot.strip!
    elsif line =~/action/
      action = line.split(/=/)[1].gsub(/\"/,"").gsub(/\;/,"")
      action.strip!
    end
  }


  # fill in the form, then submit it
  dform = page.forms.first  
  dform.rowids=rowids
  dform.annot=annot
  dform.action=action
  page2 = mech.submit(dform)

  #parse out the download link
  link = page2.body.match(/UserDownload\/\w+.txt/)
  raise "Table link not found" if (link.nil? || link == "")

  link2 = 'http://david.abcc.ncifcrf.gov/'
  link2 << "#{link}"
  puts "Table Link:"
  puts link2

  davidHash = {}
  header = []
  count = 0

  return mech.get_file(link2)
end

#------------------------------------------------------------------------------
# takes a list of refseq gene names, returns a list of functional annotations
#
# inputs: an array of gene names [TP53,MDM2,CDKN2A]
#         (optional) species name/code - default human
#         (optional) an array of annotation types to retrieve
#             see: http://david.abcc.ncifcrf.gov/content.jsp?file=DAVID_API.html
#         (optional) custom type of gene names being input (default "OFFICIAL_GENE_SYMBOL")
#
# returns a hash of geneName -> annotation type -> array of annotations
#
# example hash returned:
# {"CDKN2A" -> 
#              {"GOTERM_BP_2"  -> [GO:0005622~intracellular, 
#                                  GO:0005737~cytoplasm, 
#                                  GO:0031981~nuclear lumen]
#              }
#
#              {"KEGG_PATHWAY" -> [hsa04115:p53 signaling pathway, 
#                                  hsa05212:Pancreatic cancer, 
#                                  hsa05214:Glioma]
#              }
# }
#

def getFuncAnnotations(geneList,species="9606:Homo sapiens",customAnnos=nil,customType=nil)
  #create a DAVID url  
  url = "http://david.abcc.ncifcrf.gov/api.jsp?"
  url << "type=OFFICIAL_GENE_SYMBOL"
  url << "&ids=#{geneList.join(",")}"
  url << "&tool=annotationReport"
  
  if customAnnos.nil?
    annotations = [
                 "GOTERM_BP_3", #only GO terms from level 3 and lower
                 "GOTERM_BP_4", #not interested in matching "biological process"
                 "GOTERM_BP_5",
                 "GOTERM_CC_3",
                 "GOTERM_CC_4",
                 "GOTERM_CC_5",
                 "GOTERM_MF_3",
                 "GOTERM_MF_4",
                 "GOTERM_MF_5",
                 "KEGG_PATHWAY", #two pathway DBs
                 "BIOCARTA"]
  else
    annotations = defaultAnnos
  end


  url << "&annot=" << annotations.join(",")

  #retrieve the table
  table = retrieveTable(url)

  #parse the table
  count = 0
  header = nil
  funcHash = {}

  table.each{|line|    
    if count == 0
      header = line.chomp.split(/\t/)
    elsif line =~ /#{species}/
      data = line.chomp.split(/\t/)     

      colsDone = []

      #these three get special treatment
      pos = header.index("GENE_SYMBOL")
      name = data[pos]
      funcHash[name] = {}
      colsDone << pos      

      pos = header.index("Gene Name")
      funcHash[name]["Gene Name"] = data[pos]
      colsDone << pos

      pos = header.index("Species")
      funcHash[name]["Species"] = data[pos]
      colsDone << pos
    
      #now all the rest, comma sep lists -> arrays
      data.each_index{|i|
        unless colsDone.include?(i)
          funcHash[name][header[i]] = data[i].split(/,/).collect{|x| x.strip}
          funcHash[name][header[i]].delete_if {|x| x == "" }
        end
      }
    end
    count += 1
  }
  return funcHash
end #def


#------------------------------------------------------------------------------
# finds the functional annotations that are common between 
# all genes in the davidHash
#
# input:  a davidHash (results of getFuncAnnotations)
# output: a hash of fields containing arrays of 
#          common annotations
#
# {"GOTERM_BP_2"  -> [GO:0005622~intracellular, 
#                     GO:0005737~cytoplasm]
# }                               
#

def getCommonAnnos(davidHash)
  commonHash = {}

  fieldNames = davidHash[davidHash.keys.first].keys
  genes = davidHash.keys

  fieldNames.each{|field|
    annos = {}
    genes.each{|gene|
      if davidHash[gene].key?(field)
        davidHash[gene][field].each{|anno|
          anno = "#{anno.strip}"
          if annos.key?("#{anno}")
            annos["#{anno}"] += 1
          else
            annos["#{anno}"] = 1
          end
        }      
        annos.each{|k,v|
          #is the annotation found for every gene?
          if v == genes.length
            commonHash[field] = [] unless commonHash.key?(field)
            commonHash[field] << k
          end
        }
      end
    }
  }
  return commonHash
end

#------------------------------------------------------------------------------
# retrieves the functional annotation clustering results
# for the given gene list (uses default values, since 
# those options aren't exposed through the API )
#
# inputs: an array of gene names [TP53,MDM2,CDKN2A]
#         (optional) an array of annotation types to cluster
#             see: http://david.abcc.ncifcrf.gov/content.jsp?file=DAVID_API.html
#         (optional) custom type of gene names being input (default "OFFICIAL_GENE_SYMBOL")
#
# returns a hash of terms -> hash of info for that term
#
# example:
# {"GO:0005488~binding"  -> 
#                           {"category" ->  "GOTERM_MF_ALL",
#                            "count"    ->  "4"
#                            "PValue"   ->  "0.00000351590729492051"
#                            . . . 

def getFuncClustering(geneList,customAnnos=nil,customType=nil)

  #create a DAVID url  
  url = "http://david.abcc.ncifcrf.gov/api.jsp?"

  if customType.nil?
    url << "type=OFFICIAL_GENE_SYMBOL"
  else
    url << "type=#{customType}"
  end

  url << "&ids=#{geneList.join(",")}"
  url << "&tool=term2term"
 
  if customAnnos.nil?
    annotations = ["GOTERM_BP_3", #only GO terms from level 3 and lower
                   "GOTERM_BP_4", #not interested in matching "biological process"
                   "GOTERM_BP_5",
                   "GOTERM_CC_3",
                   "GOTERM_CC_4",
                   "GOTERM_CC_5",
                   "GOTERM_MF_3",
                   "GOTERM_MF_4",
                   "GOTERM_MF_5",
                   "KEGG_PATHWAY", #two pathway DBs
                   "BIOCARTA"]
  else
    annotations = defaultAnnos
  end

  url << "&annot=" << annotations.join(",")


  #retrieve the data from DAVID
  table = retrieveTable(url)

  #parse the table
  clusterHash = {}
  count = 0
  header = []
  table.each{|line|    
    if count == 0  #summary info, not returned for now
      summary = line.chomp.split(/\t/)

    elsif count == 1 #get header line
      header = line.chomp.split(/\t/)  

    else

      #extract the term
      data = line.chomp.split(/\t/)     
      pos = header.index("Term")
      name = data[pos]
      clusterHash[name] = {}
      
      #now add the rest of the terms
      data.each_index{|i|
        unless i == pos
          #convert to comma-sep list, if appropriate
          if data[i] =~ /,/
            data[i] = data[i].split(/,/).collect{|x| x.strip}
            #remove blank entries
            data[i] = data[i].delete_if{|x| x == ""}
          end
          
          clusterHash[name][header[i]] = data[i]
        end
      }
    end
    count += 1
  }

  return clusterHash
end #def



##------------------------------------------------------------------------------
# retrieves full gene information for the given gene list
#
# inputs: an array of gene names [TP53,MDM2,CDKN2A]
#         (optional) species name/code - default human
#         (optional) custom type of gene names being input (default "OFFICIAL_GENE_SYMBOL")
#
#
# returns a hash of genes -> hash of info for that term
#
# example:
# {"MDM2"   -> 
#               {"OMIM_DISEASE"   -> "blah blah blah"
#               {"ENTREZ_GENE_ID" -> "blah,blah,blah"
#                            . . . 
#

def getGeneReport(geneList,species="9606:Homo sapiens",customType=nil)

  #create a DAVID url  
  url = "http://david.abcc.ncifcrf.gov/api.jsp?"

  if customType.nil?
    url << "type=OFFICIAL_GENE_SYMBOL"
  else
    url << "type=#{customType}"
  end

  url << "&ids=#{geneList.join(",")}"
  url << "&tool=geneReportFull"
 

  #retrieve the data from DAVID
  table = retrieveTable(url)

  #parse the table
  #there's a more elegant way to parse this, but this'll
  #do for now (also a less mem-intensive way)
  newAnno = true
  blocks = []
  current = -1
  table.each{|line|    
    if line.chomp == ""
      newAnno = true      
    elsif newAnno #first line in block
      current += 1
      blocks[current] = []      
      blocks[current] << line.chomp.split(/\t/)
      newAnno = false
    else
      blocks[current] << line.chomp.split(/\t/)
    end
  }

  infoHash = {}
  blocks.each{|block|
    if block[0][2] =~ /#{species}/
      name = block[0][0]
      infoHash[name] = {}
      infoHash[name]["species"] = block[0][2]
      infoHash[name]["description"] = block[0][1]

      1.upto(block.length-1){|i|
        line = block[i]
        infoHash[name][line[0]] = line[1]
      }
    end
  }

  return infoHash
end
