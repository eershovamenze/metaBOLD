# metaBOLD
## Overview
This project aims to maintain a dynamically growing database of Molecular Operational Taxonomic Units (MOTUs) derived from invertebrate community metabarcoding data, based on the COI-mini barcode (Leray fragment). The concept is borrowed from BOLD's BIN (Barcode Identifier Index) principle (https://portal.boldsystems.org/bin). The core idea is to store a “reference” collection of MOTUs, each represented by a known reference sequence or consensus sequence, and assign them a global identifier index (MIN, MOTU Identifier Index). When new sequences are introduced, the system compares them against existing MOTUs. If a new sequence meets a predefined similarity threshold (e.g., 97% sequence similarity), it is assigned to that MOTU. Otherwise, it generates a new MOTU identifier (MIN) and adds it to the database.

This approach ensures that as more sequencing data are produced, the database automatically incorporates these novel sequences into a cohesive, standardized MOTU framework.

## Key Features

1. **Global data platform for storing MOTU tables**
    - At this stage, only MOTU or ASV tables of Leray COI data generated from invertebrate community samples are eligible to be added to the database (not eDNA, gut contents, or high-throughput barcoding datasets)
    - Only samples with a minimum sequencing depth of 30,000 reads are incorporated (OK?)
    - Four files are submitted:
        - File 1: XXXX.csv. Template file that includes metadata for the project. Habitat (marine, freshwater, soil, etc.); target community type (zooplankton, benthos, insects, etc.); region; researcher names/affiliations.
        - File 2: XXXX_MOTU.csv. This table contains 3 mandatory columns: MOTU_id, MOTU_sequence, and at least one sample_id column representing read counts
        - File 3: XXXX_PCR.csv. This table contains 3 mandatory columns: sample_id (corresponding to sample_id column in MOTU table), station_id, replicate, sample_type (values: "sample" or "control"). If data is pre-processed (i.e., controls are removed and replicates are pooled), then just list "1" as replicate and "sample" for sample_type for all sample_id values.
        - File 4: XXXX_stations.csv. Station metadata table (.csv format). This table includes all station metadata. Mandatory columns: station_id (corresponding to station_id from Sample), date, latitude, longitude. Optional columns: depth, sampling gear, etc.
        

2. **Data filtration**
    - Incoming data is pre-filtered to remove bacterial sequences and potential NUMTS (nuclear pseudo-genes) (how?)

3. **Dynamic MOTU assignment**:  
   - New sequences are compared against the list of existing MINs and those that are within 2(3?) similarity of an existing MIN is assigned to that MIN
   - When no match is found, a new MIN is created 
   - To limit volume, at this point only sequences that contribute a minimum of %% are processed (but all are retained in original datasets)
   - Provide stable, reproducible MOTU assignments that can be tracked across different datasets and sampling periods.


3. **Reference database maintenance**:  
   - Maintain a reference database of MOTUs, each represented by a representative sequence (or a consensus sequence)
   - Keep track of variability within MOTUs.
   - Facilitate easy updates, ensuring that as MOTU definitions evolve, historical assignments remain accessible and comparable.

4. **Associated metadata and distribution insights**:  
   - Explore geographical distributions and community compositions by querying MOTUs alongside their metadata.
   - Enable researchers to visualize how MOTUs are spread across regions and within particular communities over time.

5. **Integration with Web and SQL**:  
   - Utilize a SQL database for structured, scalable data storage.
   - Provide a web interface for searching, browsing, and visualizing MOTU distributions and community assemblages.
   - Offer an intuitive platform for researchers, lab members, and collaborators to access and update MOTU information.

## Action plan:
- Populate "initial" MOTU list
- Test filtering methods
- Test clustering methods
- Set up database and an easy way to populate it
- Set up basic web interface
- Etc.?

---

By providing a stable, scalable platform for MOTU classification and distribution analysis, this project фшьы to support researchers in tracking biodiversity shifts, understanding invertebrate community composition, and making informed conservation decisions.
