
# Get all JSON file paths in the directory
import json
import os
from glob import glob
from itertools import combinations

# Directory containing JSON files
json_dir = "."  # Change this to your directory containing the JSON files

# Get all JSON file paths in the directory

json_files = ["sigmod2022.json", "sigmod2023.json"]
json_files = glob(os.path.join(json_dir, "*.json"))
# Initialize mappings and lists
author_id_map = {}  # Maps original author IDs to integer numbers
author_name_map = {}  # Maps integer numbers to author names
author_counter = 0
coauthor_list = []

# Process each JSON file
for json_file in json_files:
    with open(json_file, "r") as file:
        data = json.load(file)

    # Extract the list of papers
    papers = data["result"]["hits"]["hit"]

    for paper in papers:
        info = paper.get("info", {})
        authors = info.get("authors", {}).get("author", [])

        # Ensure authors are in list format
        if isinstance(authors, dict):  # If only one author exists
            authors = [authors]

        # Assign integer IDs to authors
        author_numbers = []
        for author in authors:
            pid = author["@pid"]
            name = author["text"]
            if pid not in author_id_map:
                author_id_map[pid] = author_counter
                author_name_map[author_counter] = name
                author_counter += 1
            author_numbers.append(str(author_id_map[pid]))

        # Store co-author integer IDs (comma-separated) if there are more than one
        if len(author_numbers) > 1:
            coauthor_list.append(author_numbers)

# Generate edge list from co-author groups
edges = []
for coauthors in coauthor_list:
    for edge in combinations(coauthors, 2):
        edges.append(" ".join(edge))

# Save edge list to a CSV file
with open("edge_list.csv", "w") as f:
    for edge in edges:
        f.write(edge + "\n")

# Save author mapping to a CSV file
with open("author_mapping.csv", "w") as f:
    f.write("Integer ID, Author Name\n")  # CSV header
    for int_id in sorted(author_name_map.keys()):
        f.write(f"{int_id}, {author_name_map[int_id]}\n")

print(f"Processed {len(json_files)} files.")
print(f"Generated {len(edges)} edges and mapped {len(author_name_map)} authors.")
