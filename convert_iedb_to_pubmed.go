package iedb_to_pubmed

import (
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"regexp"
	"sort"
	"strings"
)

var (
	REintrafield_delim = regexp.MustCompile(`,\s+`)

	current_epitope_id string = "-1"
	pss                map[string]int
	pmids              map[string]int
	smiles             string
)

func splitify(v string, store map[string]int) {
	if v == "" {
		return
	}

	if REintrafield_delim.MatchString(v) {
		for _, s := range strings.Split(v, ", ") {
			if s != "null" {
				store[s] = 1
			}
		}
	} else if v != "" && v != "null" {
		store[v] = 1
	}
}

func pubmedify(epitope_id string, pss map[string]int, smiles string, pmids map[string]int, outstream *os.File) {
	fmt.Fprint(outstream, strings.Join([]string{
		"Epitope ID:" + epitope_id,
		newlineify(pss),
		"www.iedb.org",
		"http://www.iedb.org/epId/" + epitope_id,
		smiles,
		newlineify(pmids),
	}, ",")+"\r\n")
}

func newlineify(collection map[string]int) string {
	strs := make([]string, 0)
	for k, _ := range collection {
		strs = append(strs, k)
	}
	sort.Strings(strs)
	return `"` + strings.Join(strs, "\n") + `"`
}

func lineHandler(row []string, outstream *os.File) {

	accession := row[1]
	// reject rows without CHEBI: in accession
	if !strings.HasPrefix(accession, "CHEBI:") {
		return
	}

	epitope_id := row[0]
	if current_epitope_id != epitope_id {
		if current_epitope_id != "-1" {
			pubmedify(current_epitope_id, pss, smiles, pmids, outstream)
		}
		current_epitope_id = epitope_id
		pss = make(map[string]int)
		pmids = make(map[string]int)
		smiles = ""
	}

	aliases := row[2]
	synonyms := row[3]
	pubchem_pubmed_id := row[5]

	smiles = row[4]

	if accession != "" {
		pss[accession] = 1
	}
	splitify(aliases, pss)
	splitify(synonyms, pss)
	if pubchem_pubmed_id != "" {
		pmids[pubchem_pubmed_id] = 1
	}
}

func process_file(path string) {
	file, err := os.Open(path)
	if err != nil {
		fmt.Printf("error opening file: %v\n", err)
		os.Exit(1)
	}

	data := csv.NewReader(file)
	fmt.Fprint(os.Stdout, "PUBCHEM_EXT_DATASOURCE_REGID,PUBCHEM_SUBSTANCE_SYNONYM,PUBCHEM_EXT_DATASOURCE_URL,PUBCHEM_EXT_SUBSTANCE_URL,PUBCHEM_EXT_DATASOURCE_SMILES,PUBCHEM_PUBMED_ID\r\n")
	for {
		row, rerr := data.Read()
		if rerr == io.EOF {
			pubmedify(current_epitope_id, pss, smiles, pmids, os.Stdout)
			break
		}
		if rerr != nil {
			if rerr != io.EOF {
				fmt.Printf("error reading file: %v\n", rerr)
				os.Exit(1)
			}
		}
		if row == nil || row[0] == "epitope_id" {
			continue
		}
		lineHandler(row, os.Stdout)
	}
}
