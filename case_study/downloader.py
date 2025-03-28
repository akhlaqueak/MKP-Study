import requests
import sys
import os

conf=sys.argv[1]
os.makedirs(conf, exist_ok=True)
# Define the range of years and base URL
start_year = 2004
end_year = 2022
base_url = f"https://dblp.org/search/publ/api?q=toc%3Adb/conf/{conf}/{conf}.bht%3A&h=1000&format=json"

# Loop through the years and download JSON files
for year in range(start_year, end_year + 1):
    # url = f"https://dblp.org/search/publ/api?q=toc%3Adb/conf/icde/icde{year}.bht%3A&h=1000&format=json"
    url = f"https://dblp.org/search/publ/api?q=toc%3Adb/journals/vldb/vldb{year}.bht%3A&h=1000&format=json"
    url = f"https://dblp.org/search/publ/api?q=toc%3Adb/conf/sigmod/sigmod{year}.bht%3A&h=1000&format=json"
    # url = base_url.format(year)
    print(url)
    response = requests.get(url)
    
    if response.status_code == 200:
        filename = f"{conf}/{conf}_{year}.json"
        with open(filename, "w", encoding="utf-8") as f:
            f.write(response.text)
        print(f"Downloaded: {filename}")
    else:
        print(f"Failed to download ICDE {year} data. HTTP Status: {response.status_code}")
