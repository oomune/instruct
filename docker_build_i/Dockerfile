# setting base image
FROM continuumio/anaconda:latest

ADD tsv2rdf_instruct.py /tmp
ADD templ_instruct.ttl /tmp
ADD templ_instruct.ttl.prefix /tmp
ADD templ_instruct.ttl.evi /tmp
ENTRYPOINT ["python", "/tmp/tsv2rdf_instruct.py"]
CMD ["--c","/mnt/tsv2rdf_instruct.json"]
