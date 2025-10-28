import tables
from pathlib import Path
from dataclasses import dataclass
import typing as tp

@dataclass
class TableAttribute:
    name: str
    predicate: tp.Callable = lambda x: True
    transform: tp.Callable = lambda x: x

def get_table_attributes(cls):
    attrs = {}
    for field_name, annotated_type in tp.get_type_hints(cls, include_extras=True).items():
        if tp.get_origin(annotated_type) is tp.Annotated:
            _, table_attr = tp.get_args(annotated_type)
            attrs[field_name] = table_attr
    return attrs

class H5Reader:
    def __init__(self, h5_file: Path) -> None:
        if not h5_file.exists():
            raise FileNotFoundError(f"The file {h5_file} does not exist.")
        self.h5_file = h5_file
        self.handle: tables.File | None = None 

    def open(self) -> None:
        self.handle = tables.open_file(str(self.h5_file), mode='r')
    
    def close(self) -> None:
        if self.handle:
            self.handle.close()
            self.handle = None
    
    def __enter__(self) -> "H5Reader":
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def _check(self, obj: tp.Any) -> dict[str, TableAttribute]:
        attrs = get_table_attributes(obj)
        if not attrs:
            raise ValueError(f"No TableAttributes found in {obj.__name__}.")
        if "_table" not in obj.__dict__:
            raise ValueError(f'Missing required TableAttribute "_table" in {obj.__name__}.')
        return obj._table, attrs
    
    def _apply_override(self, attrs, override):
        if not override:
            return
        
        for field_name, attr in attrs.items():
            if field_name in override:
                if "transform" in override[field_name]:
                    attr.transform = override[field_name]["transform"]
                if "predicate" in override[field_name]["predicate"]:
                    attr.predicate = override[field_name]["predicate"]

    def fetch(self, obj: tp.Any, override: dict = {}) -> tp.Generator[tp.Any, None, None]:
        if self.handle is None:
            raise RuntimeError("H5 file is not opened. Use 'with' statement or call open() method.")
        
        table_name, attrs = self._check(obj)
        h5_table = self.handle.get_node(table_name)
        self._apply_override(attrs, override)

        for row in h5_table.iterrows():
            args = {}
            for field_name, attr in attrs.items():
                value = row[attr.name]
                if not attr.predicate(value):
                    break
                else:
                    args[field_name] = attr.transform(value)
            else:
                yield obj(**args)

    def fetch_raw(self, obj: tp.Any) -> tp.Generator[tp.Any, None, None]:
        if self.handle is None:
            raise RuntimeError("H5 file is not opened. Use 'with' statement or call open() method.")
        
        table_name, attrs = self._check(obj)
        h5_table = self.handle.get_node(table_name)
        
        for row in h5_table.iterrows():
            obj(**{field_name: row[attr.name] for field_name, attr in attrs.items()})
    
    def fetch_rows(self, obj: tp.Any) -> tp.Generator[tp.Any, None, None]:
        if self.handle is None:
            raise RuntimeError("H5 file is not opened. Use 'with' statement or call open() method.")
        table_name, _ = self._check(obj)
        yield from self.handle.get_node(table_name).iterrows()
        
