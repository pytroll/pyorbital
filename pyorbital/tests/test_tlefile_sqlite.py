"""Test TLE SQLiteTLE."""

import datetime as dt
import sqlite3
import time
from pathlib import Path

import pytest

from pyorbital.tests.test_tlefile import LINE1, LINE1_2, LINE2, LINE2_2
from pyorbital.tlefile import (
    ISO_TIME_FORMAT,
    PLATFORM_NAMES_TABLE,
    SATID_TABLE,
    ChecksumError,
    SQLiteTLE,
    Tle,
    _utcnow,
    table_exists,
)


@pytest.fixture
def platforms():
    """Provide a default ISS platform mapping."""
    return {25544: "ISS"}


@pytest.fixture
def writer_config(tmp_path):
    """Provide writer configuration for TLE output."""
    return {
        "output_dir": str(tmp_path / "tle_dir"),
        "filename_pattern": "tle_%Y%m%d_%H%M%S.%f.txt",
        "write_name": True,
        "write_always": False,
    }


@pytest.fixture
def db(tmp_path, platforms, writer_config):
    """Create a SQLiteTLE instance backed by a temporary database."""
    db_path = tmp_path / "tle.db"
    db = SQLiteTLE(str(db_path), platforms, writer_config)
    yield db
    db.close()


@pytest.fixture
def tle():
    """Provide a valid ISS TLE instance."""
    return Tle("ISS", line1=LINE1, line2=LINE2)


def test_init(db):
    """Verify that the platform_names table is created with the correct schema."""
    columns = [col.strip() for col in PLATFORM_NAMES_TABLE.strip("()").split(",")]
    num_columns = len(columns)

    assert table_exists(db.db, "platform_names")

    res = db.db.execute("SELECT * FROM platform_names")
    names = [desc[0] for desc in res.description]

    assert len(names) == num_columns
    for col in columns:
        assert col.split(" ")[0] in names


def test_update_db(db, tle, platforms):
    """Verify update_db creates tables, inserts rows, and ignores duplicates."""
    satid = list(platforms.keys())[0]

    # Case 1: platform not configured
    db.platforms = {}
    db.update_db(tle, "foo")
    assert not table_exists(db.db, satid)
    assert not db.updated

    # Case 2: platform configured
    db.platforms = platforms
    db.update_db(tle, "foo")
    assert table_exists(db.db, satid)
    assert db.updated

    # Check columns
    columns = [col.strip() for col in SATID_TABLE.replace("'{}' (", "").strip(")").split(",")]
    res = db.db.execute(f"SELECT * FROM '{satid:d}'")
    names = [desc[0] for desc in res.description]
    for col in columns:
        assert col.split(" ")[0] in names

    # Check data
    data = res.fetchall()
    assert len(data) == 1
    assert data[0][0] == "2008-09-20T12:25:40.104192"
    assert data[0][1] == "\n".join((LINE1, LINE2))

    date_added = dt.datetime.strptime(data[0][2], ISO_TIME_FORMAT)
    assert (_utcnow() - date_added).total_seconds() < 1.0
    assert data[0][3] == "foo"

    # Duplicate epoch → no new row
    db.update_db(tle, "bar")
    res = db.db.execute(f"SELECT * FROM '{satid:d}'")
    data2 = res.fetchall()
    assert len(data2) == 1
    assert data2[0][2] == data[0][2]
    assert data2[0][3] == "foo"


def test_update_db_newer_epoch(db, tle):
    """Verify update_db inserts a new row when the epoch is newer."""
    db.update_db(tle, "foo")

    tle2 = Tle("ISS", line1=LINE1, line2=LINE2)

    class Dummy:
        def __init__(self, dt):
            self._dt = dt

        def item(self):
            return self._dt

    tle2.epoch = Dummy(tle.epoch + dt.timedelta(seconds=10))

    db.update_db(tle2, "bar")

    rows = db.db.execute("SELECT * FROM '25544' ORDER BY epoch").fetchall()
    assert len(rows) == 2


def test_update_db_invalid_tle(db):
    """Verify that constructing an invalid TLE raises an exception."""
    BAD_LINE1 = "1 25544U 98067A   00000.00000000  .00000000  00000-0  00000-0 0  0000"
    BAD_LINE2 = "2 25544  00.0000 000.0000 0000000 000.0000 000.0000 00.00000000000000"

    with pytest.raises(ChecksumError):
        Tle("ISS", line1=BAD_LINE1, line2=BAD_LINE2)


def test_write_tle_txt(db, tle, writer_config):
    """Verify write_tle_txt writes files according to update and config rules."""
    out_dir = Path(writer_config["output_dir"])

    db.update_db(tle, "foo")

    # Case 1: updated=False → no file
    db.updated = False
    db.write_tle_txt()
    assert not out_dir.exists()

    # Case 2: updated=True → file written
    db.updated = True
    db.write_tle_txt()
    assert out_dir.exists()

    files = list(out_dir.glob("tle_*txt"))
    assert len(files) == 1
    assert "%" not in files[0].name

    content = files[0].read_text().splitlines()
    assert content == ["ISS", LINE1, LINE2]

    # Case 3: updated=False again → no new file
    db.updated = False
    db.write_tle_txt()
    assert len(list(out_dir.glob("tle_*txt"))) == 1

    # Case 4: write_always=True, write_name=False → new file
    db.writer_config["write_always"] = True
    db.writer_config["write_name"] = False

    time.sleep(1.5)
    db.write_tle_txt()

    files = sorted(out_dir.glob("tle_*txt"))
    assert len(files) == 2

    content = files[1].read_text().splitlines()
    assert content == [LINE1, LINE2]


def test_write_tle_txt_missing_table(tmp_path, writer_config):
    """Verify write_tle_txt fails when a satellite table is missing."""
    db_path = tmp_path / "tle.db"
    db = SQLiteTLE(str(db_path), {99999: "MISSING"}, writer_config)
    db.updated = True

    with pytest.raises(sqlite3.OperationalError):
        db.write_tle_txt()

    db.close()


def test_write_tle_txt_empty_table(tmp_path, writer_config):
    """Verify write_tle_txt fails when a satellite table exists but has no rows."""
    db_path = tmp_path / "tle.db"
    db = SQLiteTLE(str(db_path), {25544: "ISS"}, writer_config)

    db.db.execute(f"CREATE TABLE {SATID_TABLE.format(25544)}")

    db.updated = True
    with pytest.raises(TypeError):
        db.write_tle_txt()

    db.close()


def test_nested_filename_pattern(tmp_path):
    """Verify nested strftime patterns in output_dir create the expected directories."""
    writer_config = {
        "output_dir": str(tmp_path / "out/%Y/%m"),
        "filename_pattern": "tle_%H%M%S.txt",
        "write_name": True,
        "write_always": True,
    }

    db_path = tmp_path / "tle.db"
    db = SQLiteTLE(str(db_path), {25544: "ISS"}, writer_config)
    tle = Tle("ISS", line1=LINE1, line2=LINE2)

    db.update_db(tle, "foo")
    db.write_tle_txt()

    expanded_dir = Path(_utcnow().strftime(writer_config["output_dir"]))
    assert expanded_dir.exists()

    files = list(expanded_dir.glob("tle_*txt"))
    assert len(files) == 1

    db.close()


def test_write_multiple_satellites(tmp_path):
    """Verify write_tle_txt writes TLEs for multiple satellites in order."""
    platforms = {25544: "ISS", 38771: "SAT38771"}

    writer_config = {
        "output_dir": str(tmp_path / "tle_dir"),
        "filename_pattern": "tle_%Y%m%d_%H%M%S.%f.txt",
        "write_name": True,
        "write_always": True,
    }

    db_path = tmp_path / "tle.db"
    db = SQLiteTLE(str(db_path), platforms, writer_config)

    tle1 = Tle("ISS", line1=LINE1, line2=LINE2)
    tle2 = Tle("SAT38771", line1=LINE1_2, line2=LINE2_2)

    db.update_db(tle1, "foo")
    db.update_db(tle2, "bar")

    db.write_tle_txt()

    out_dir = Path(writer_config["output_dir"])
    files = list(out_dir.glob("tle_*txt"))
    assert len(files) == 1

    content = files[0].read_text().splitlines()
    assert content == [
        "ISS",
        LINE1,
        LINE2,
        "SAT38771",
        LINE1_2,
        LINE2_2,
    ]

    db.close()


def test_close_closes_connection(tmp_path, platforms, writer_config):
    """Verify that close() prevents further database operations."""
    db_path = tmp_path / "tle.db"
    db = SQLiteTLE(str(db_path), platforms, writer_config)

    assert table_exists(db.db, "platform_names")

    db.close()

    import sqlite3

    with pytest.raises(sqlite3.ProgrammingError):
        db.db.execute("SELECT 1")


SATELLITES = [
    ({25544: "ISS"}, "ISS", LINE1, LINE2),
    ({38771: "SAT38771"}, "SAT38771", LINE1_2, LINE2_2),
]


@pytest.mark.parametrize(("platforms", "expected_name", "line1", "line2"), SATELLITES)
def test_update_db_platforms(tmp_path, writer_config, platforms, expected_name, line1, line2):
    """Verify update_db works for multiple satellites with matching TLEs."""
    db_path = tmp_path / "tle.db"
    db = SQLiteTLE(str(db_path), platforms, writer_config)

    tle = Tle(expected_name, line1=line1, line2=line2)
    satid = list(platforms.keys())[0]

    db.update_db(tle, "foo")
    assert table_exists(db.db, satid)

    res = db.db.execute(f"SELECT * FROM '{satid:d}'")
    data = res.fetchall()
    assert len(data) == 1
    assert data[0][1] == "\n".join((line1, line2))

    db.close()
