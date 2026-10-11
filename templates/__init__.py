import os

from fastapi.templating import Jinja2Templates

from util import relative_path

templates = Jinja2Templates(directory=relative_path("templates"))

# added to the urls of our own javascript and css, so a browser doesn't keep using the old ones after an update
templates.env.globals["static_version"] = int(max(
    os.path.getmtime(relative_path(os.path.join("static", name))) for name in ("molclass.js", "style.css")
))
