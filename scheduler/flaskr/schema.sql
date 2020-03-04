DROP TABLE IF EXISTS acc;
DROP TABLE IF EXISTS fastq;
DROP TABLE IF EXISTS acc_state_defn;
DROP TABLE IF EXISTS fastq_state_defn;

CREATE TABLE acc_state_defn(
    acc_state_id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT
);

CREATE TABLE fastq_state_defn(
    fastq_state_id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT
);

CREATE TABLE acc(
    acc_id INTEGER PRIMARY KEY AUTOINCREMENT,
    acc_state_id INTEGER,

    acc TEXT,
    sra_url TEXT,
    split_cmd TEXT,
    merge_cmd TEXT,
    FOREIGN KEY(acc_state_id) REFERENCES acc_state_defn(acc_state_id)
);

CREATE TABLE fastq(
    fastq_id INTEGER PRIMARY KEY AUTOINCREMENT,
    fastq_state_id INTEGER,
    acc_id INTEGER,
    n INTEGER,
    align_cmd TEXT,
    paired BOOLEAN,
    FOREIGN KEY(acc_id) REFERENCES acc(acc_id)
    FOREIGN KEY(fastq_state_id) REFERENCES fastq_state_defn(fastq_state_id)
);
