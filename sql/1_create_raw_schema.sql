------------------------------------------------------------
-- Data source directory list
------------------------------------------------------------
CREATE TABLE IF NOT EXISTS data_source (
    exp_code TEXT PRIMARY KEY         -- directory / file prefix
);

------------------------------------------------------------
-- MOTU table: required motu_id + sequence
------------------------------------------------------------
CREATE TABLE IF NOT EXISTS raw.motu (
    exp_code   TEXT REFERENCES data_source,
    motu_id    TEXT,
    sequence   TEXT NOT NULL,
    PRIMARY KEY (exp_code, motu_id)
);

------------------------------------------------------------
-- Samples table: sample_id belongs to a data_source
------------------------------------------------------------
CREATE TABLE IF NOT EXISTS raw.samples (
    exp_code    TEXT REFERENCES data_source,
    sample_id   TEXT,                 -- column names in reads matrix
    type        TEXT,                 -- 'sample' | 'blank'
    sample_name TEXT,                 -- should match station table
    PRIMARY KEY (exp_code, sample_id)
);

------------------------------------------------------------
-- Stations table: coordinates keyed by sample_name
------------------------------------------------------------
CREATE TABLE IF NOT EXISTS raw.stations (
    exp_code    TEXT REFERENCES data_source,
    sample_name TEXT,
    latitude    NUMERIC,
    longitude   NUMERIC,
    extra       JSONB,                -- optional, variable columns
    PRIMARY KEY (exp_code, sample_name)
);

------------------------------------------------------------
-- Reads table (long form) : one row per MOTU × sample
------------------------------------------------------------
CREATE TABLE IF NOT EXISTS raw.reads (
    exp_code   TEXT,
    motu_id    TEXT,
    sample_id  TEXT,
    read_count INTEGER NOT NULL CHECK (read_count >= 0),
    PRIMARY KEY (exp_code, motu_id, sample_id),
    FOREIGN KEY (exp_code, motu_id)   REFERENCES raw.motu    (exp_code, motu_id),
    FOREIGN KEY (exp_code, sample_id) REFERENCES raw.samples (exp_code, sample_id)
);

CREATE INDEX IF NOT EXISTS raw_reads_motu_idx  ON raw.reads (motu_id);
CREATE INDEX IF NOT EXISTS raw_reads_sample_idx ON raw.reads (sample_id);
