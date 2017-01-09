--
-- PostgreSQL database dump
--

SET statement_timeout = 0;
SET lock_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SET check_function_bodies = false;
SET client_min_messages = warning;

--
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


--
-- Name: plv8; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plv8 WITH SCHEMA pg_catalog;


--
-- Name: EXTENSION plv8; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plv8 IS 'PL/JavaScript (v8) trusted procedural language';


SET search_path = public, pg_catalog;

SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: idmapping; Type: TABLE; Schema: public; Owner: agoncear; Tablespace: 
--

CREATE TABLE idmapping (
    uniprot_ac character varying NOT NULL,
    uniprot_id character varying,
    gene_id integer,
    gi integer,
    refseq character varying,
    pdb character varying,
    go character varying,
    ensembl character varying
);


ALTER TABLE idmapping OWNER TO agoncear;

--
-- Name: pairs; Type: TABLE; Schema: public; Owner: agoncear; Tablespace: 
--

CREATE TABLE pairs (
    pair_id integer NOT NULL,
    tax_id integer NOT NULL,
    rev integer NOT NULL,
    id_a character varying NOT NULL,
    id_b character varying NOT NULL,
    uniprot_ac_a character varying NOT NULL,
    uniprot_ac_b character varying NOT NULL,
    tpl character varying,
    query_type character varying,
    template_type character varying,
    score_model_z double precision,
    score_model_minus_avg double precision,
    score3 double precision,
    score4 double precision,
    score5 double precision,
    score6 double precision,
    identical_a integer,
    identical_b integer,
    positive_a integer,
    positive_b integer,
    aln_len_a integer,
    aln_len_b integer,
    bs_len_a integer,
    bs_len_b integer,
    bs_covered_a integer,
    bs_covered_b integer,
    bs_aligned_a integer,
    bs_aligned_b integer,
    bs_identical_a integer,
    bs_identical_b integer,
    bs_positive_a integer,
    bs_positive_b integer,
    bs_contacts_a integer,
    bs_contacts_b integer,
    bs_blosum_a double precision,
    bs_blosum_b double precision,
    bs_score1_a double precision,
    bs_score1_b double precision,
    site_a character varying,
    site_b character varying
);


ALTER TABLE pairs OWNER TO agoncear;

--
-- Name: pairs_pair_id_seq; Type: SEQUENCE; Schema: public; Owner: agoncear
--

CREATE SEQUENCE pairs_pair_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE pairs_pair_id_seq OWNER TO agoncear;

--
-- Name: pairs_pair_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: agoncear
--

ALTER SEQUENCE pairs_pair_id_seq OWNED BY pairs.pair_id;


--
-- Name: proteins; Type: TABLE; Schema: public; Owner: agoncear; Tablespace: 
--

CREATE TABLE proteins (
    uniprot_ac character varying NOT NULL,
    uniprot_accessions character varying,
    uniprot_id character varying,
    tax_id integer NOT NULL,
    name character varying,
    sequence character varying,
    seqlen integer
);


ALTER TABLE proteins OWNER TO agoncear;

--
-- Name: species; Type: TABLE; Schema: public; Owner: agoncear; Tablespace: 
--

CREATE TABLE species (
    tax_id integer NOT NULL,
    name character varying
);


ALTER TABLE species OWNER TO agoncear;

--
-- Name: species_tax_id_seq; Type: SEQUENCE; Schema: public; Owner: agoncear
--

CREATE SEQUENCE species_tax_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE species_tax_id_seq OWNER TO agoncear;

--
-- Name: species_tax_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: agoncear
--

ALTER SEQUENCE species_tax_id_seq OWNED BY species.tax_id;


--
-- Name: pair_id; Type: DEFAULT; Schema: public; Owner: agoncear
--

ALTER TABLE ONLY pairs ALTER COLUMN pair_id SET DEFAULT nextval('pairs_pair_id_seq'::regclass);


--
-- Name: tax_id; Type: DEFAULT; Schema: public; Owner: agoncear
--

ALTER TABLE ONLY species ALTER COLUMN tax_id SET DEFAULT nextval('species_tax_id_seq'::regclass);


--
-- Name: idmapping_pkey; Type: CONSTRAINT; Schema: public; Owner: agoncear; Tablespace: 
--

ALTER TABLE ONLY idmapping
    ADD CONSTRAINT idmapping_pkey PRIMARY KEY (uniprot_ac);


--
-- Name: pairs_pkey; Type: CONSTRAINT; Schema: public; Owner: agoncear; Tablespace: 
--

ALTER TABLE ONLY pairs
    ADD CONSTRAINT pairs_pkey PRIMARY KEY (pair_id);


--
-- Name: proteins_pkey; Type: CONSTRAINT; Schema: public; Owner: agoncear; Tablespace: 
--

ALTER TABLE ONLY proteins
    ADD CONSTRAINT proteins_pkey PRIMARY KEY (uniprot_ac);


--
-- Name: species_pkey; Type: CONSTRAINT; Schema: public; Owner: agoncear; Tablespace: 
--

ALTER TABLE ONLY species
    ADD CONSTRAINT species_pkey PRIMARY KEY (tax_id);


--
-- Name: pairs_tax_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: agoncear
--

ALTER TABLE ONLY pairs
    ADD CONSTRAINT pairs_tax_id_fkey FOREIGN KEY (tax_id) REFERENCES species(tax_id);


--
-- Name: proteins_tax_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: agoncear
--

ALTER TABLE ONLY proteins
    ADD CONSTRAINT proteins_tax_id_fkey FOREIGN KEY (tax_id) REFERENCES species(tax_id);


--
-- Name: public; Type: ACL; Schema: -; Owner: agoncear
--

REVOKE ALL ON SCHEMA public FROM PUBLIC;
REVOKE ALL ON SCHEMA public FROM agoncear;
GRANT ALL ON SCHEMA public TO agoncear;
GRANT ALL ON SCHEMA public TO PUBLIC;


--
-- PostgreSQL database dump complete
--

