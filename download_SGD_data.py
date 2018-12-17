# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The following two lines will be needed in every python script:
from intermine.webservice import Service
service = Service("https://yeastmine.yeastgenome.org:443/yeastmine/service", token="YOUR-API-KEY")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "primaryIdentifier", "secondaryIdentifier", "organism.shortName", "symbol",
    "name"
)

# You can edit the constraint values below
query.add_constraint("Gene", "IN", "Gene list for S. cerevisiae 17 Dec 2018 11.4", code="A")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A")

for row in query.rows():
    print row["primaryIdentifier"], row["secondaryIdentifier"], row["organism.shortName"], \
        row["symbol"], row["name"]
