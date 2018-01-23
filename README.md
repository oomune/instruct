# INstruct

RDF Converter for INstruct

## Downloading data
Download data from the following web site.

    http://instruct.yulab.org/downloads.html

Store the source file(tsv) to be converted in the following directory.

    data_i

## Cinfiguration
Edit the following definition file.

    tsv2rdf_instruct.json


```
{
	"organism":{
		"sapiens":{
			"all":"sapiens.sin",
			"annotations":"direct_homology.txt"
		},
		"elegans":{
			"all":"elegans.sin"
		},
		"melanogaster":{
			"all":"melanogaster.sin"
		},
		"musculus":{
			"all":"musculus.sin"
		},
		"pombe":{
			"all":"pombe.sin"
		},
		"thaliana":{
			"all":"thaliana.sin"
		}
	},
	"template":{
		"prefix":"templ_instruct.ttl.prefix",
		"body":"templ_instruct.ttl",
		"evidence":"templ_instruct.ttl.evi"
	},
	"output_file":"instruct_sample.ttl",
	"data_path":"data_i"

}
```


## Usage

    $ ./tsv2rdf_instruct.sh
