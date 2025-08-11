from IPython.display import HTML, display

def display_filtered_table(collections):
    headers = ["Collection Short Name", "CMR Virtual Directory Page"]
    rows = []

    for concept_id, short_name in sorted(collections.items(), key=lambda x: x[1]):
        url = f"https://cmr.earthdata.nasa.gov/virtual-directory/collections/{concept_id}/temporal"
        clickable_url = f'<a href="{url}" target="_blank">{url}</a>'
        rows.append([short_name, clickable_url])

    table_rows = ""
    for row in rows:
        table_rows += "<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>"

    html_content = f"""
    <div style="margin-bottom: 10px;">
        <input type="text" id="tableSearch" placeholder="Search collections..." 
               style="width: 100%; padding: 8px; font-size: 16px;" onkeyup="filterTable()">
    </div>

    <div style="height: 400px; overflow: auto; border: 1px solid #ccc; padding: 10px; font-size: 16px;">
        <table id="collectionsTable" style="width: 100%; border-collapse: collapse;">
            <thead>
                <tr>
                    <th style="text-align: left;">{headers[0]}</th>
                    <th style="text-align: left;">{headers[1]}</th>
                </tr>
            </thead>
            <tbody>
                {table_rows}
            </tbody>
        </table>
    </div>

    <script>
    function filterTable() {{
        var input = document.getElementById("tableSearch");
        var filter = input.value.toLowerCase();
        var table = document.getElementById("collectionsTable");
        var rows = table.getElementsByTagName("tr");

        for (var i = 1; i < rows.length; i++) {{
            var row = rows[i];
            var text = row.textContent.toLowerCase();
            if (text.includes(filter)) {{
                row.style.display = "";
            }} else {{
                row.style.display = "none";
            }}
        }}
    }}
    </script>
    """

    display(HTML(html_content))
