Systematic Literature Review
The literature review was implemented in R (version 3.6.1, [97,98]) using the package easyPubMed (version 2.13, [99]) and consists of 2 steps. First an automated pubmed search was carried on April 13th, 2020 out using the query string: 

"(NGS OR \"next generation sequencing\" OR RNA-Seq OR \"mRNA sequencing\" OR \"whole exome sequencing\" OR \"whole-exome sequencing\" OR \"high throughput sequencing\" OR \"high-throughput sequencing\" OR \"DNA sequencing\" OR \"RNA sequencing\" OR WXS OR WGS OR \"whole-genome sequencing\" OR \"whole genome sequencing\") AND (rheumatology OR \"rheumatologic disease\" OR \"rheumatologic disease\" ))". 

This search resulted in 814 entries. The keywords of each returning dataset were intersected with official disease names extracted from ICD-11 [97] in order to filter out keywords that are not disease names. The remaining 253 keywords were then manually inspected to find rheumatic diseases. This approach identified the following diseases: Sacroiliitis, Deficiency of adenosine deaminase 2, chilblain lupus, giant cell arteritis, neuropsychiatric lupus, lupus nephritis, Immunoglobulin A nephropathy, sarcoidosis, rheumatoid arthritis, systemic sclerosis, systemic lupus erythematosus, juvenile idiopathic arthritis, myositis, osteoarthritis, spondyloarthritis, dermatomyositis, familial mediterranean fever, gout, uveitis, vasculitis, granulomatosis with polyangiitis, psoriatic arthritis, sjögren syndrome and synovitis.
 
In a second step more specific pubmed search was carried out using the disease names identified in the first step: 

"(NGS OR \"next generation sequencing\" OR RNA-Seq OR \"mRNA sequencing\" OR \"whole exome sequencing\" OR \"whole-exome sequencing\" OR \"high throughput sequencing\" OR \"high-throughput sequencing\" OR WXS OR WGS OR \"whole-genome sequencing\" OR \"whole genome sequencing\") AND (\"Deficiency of adenosine deaminase 2\" OR \"chilblain lupus\" OR \"giant cell arteritis\" OR \"neuropsychiatric lupus\" OR \"lupus nephritis\" OR \"Immunoglobulin A nephropathy\" OR \"sarcoidosis\" OR \"rheumatoid arthritis\" OR \"systemic sclerosis\" OR \"systemic lupus erythematosus\" OR \"juvenile idiopathic arthritis\" OR \"myositis\" OR \"osteoarthritis\" OR \"spondyloarthritis\" OR \"dermatomyositis\" OR \"familial mediterranean fever\" OR \"gout\" OR \"uveitis\" OR \"vasculitis\" OR \"granulomatosis with polyangiitis\" OR \"psoriatic arthritis\" OR \"sjögren syndrome\" OR \"synovitis\")"
 
This search was carried on April 13th, 2020 and resulted in 613 pubmed hits, which were (if possible) annotated regarding disease name, pubmed ID, assay, journal, year of publication by automatic screening the title and abstract. Reviews (i.e. publications which have “Review” in metadata) and commentaries were excluded and missing information was added manually by manual inspection of the publication . After manually deleting entries (e.g. because they do not report HTS or rheumatic diseases was not the focus of the work) 402 well annotated entries remained.