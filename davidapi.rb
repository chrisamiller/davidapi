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



# takes a list of refseq gene names, returns a list of functional annotations
# this is not a completely generic method - could be extended to allow selection
# of different annotation types, inputs (such as affy probe ids, etc). For now,
# this is it, though.

# takes an array of gene names [TP53,MDM2,CDKN2A]
# returns a hash of geneName -> annotation type -> array of annotations
#
# example:
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

def getFuncAnnotations(geneList,species="9606:Homo sapiens")

  #create a DAVID url
  
  url = "http://david.abcc.ncifcrf.gov/api.jsp?"
  url << "type=OFFICIAL_GENE_SYMBOL"
  url << "&ids=#{geneList.join(",")}"
  url << "&tool=annotationReport"
  
  annotations = [
                 "GOTERM_BP_2",               
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

  url << "&annot=" << annotations.join(",")


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
  #keep only lines that match the species                          
  mech.get_file(link2).each{|line|
    if count == 0
      header = line.chomp.split(/\t/)
    elsif line =~ /#{species}/
      data = line.chomp.split(/\t/)     

      colsDone = []

      #these three get special treatment
      pos = header.index("GENE_SYMBOL")
      name = data[pos]
      davidHash[name] = {}
      colsDone << pos      

      pos = header.index("Gene Name")
      davidHash[name]["Gene Name"] = data[pos]
      colsDone << pos

      pos = header.index("Species")
      davidHash[name]["Species"] = data[pos]
      colsDone << pos
    
      #now all the rest, comma sep lists -> arrays
      data.each_index{|i|
        unless colsDone.include?(i)
          davidHash[name][header[i]] = data[i].split(/,/).collect{|x| x.strip}
          davidHash[name][header[i]].delete_if {|x| x == "" }
        end
      }
    end
    count += 1
  }
  return davidHash
end #def



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
    puts field
    annos = {}
    genes.each{|gene|
      davidHash[gene][field].each{|anno|
        anno = "#{anno.strip}"
        puts "  testing #{anno}"
        if annos.key?("#{anno}")
          puts "found"
          annos["#{anno}"] += 1
        else
          annos["#{anno}"] = 1
        end
      }      
      annos.each{|k,v|
        #is the annotation found for every gene?
        puts "#{k} -> #{v}"
        if v == genes.length
          commonHash[field] = [] unless commonHash.key?(field)
          commonHash[field] << k
        end
      }
    }
  }
  return commonHash
end
