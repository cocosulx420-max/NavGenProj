--!strict
-- NavGen.Boundary — outlines from the per-part LOCAL GRIDS.
--
-- WHAT CHANGED, AND WHY IT IS THE WHOLE DESIGN. The previous version of this
-- module read the world-aligned surfel raster from `Floor`. On that raster a
-- rotated part's rim is a jagged run of axis-aligned cells, so the module spent
-- almost all of its length undoing that: a Euclidean distance transform, a
-- greedy total-least-squares line fit to recover the angle the raster had
-- destroyed, corner recovery from where the fit failed, then miter and bevel
-- cleanup. Every one of those stages was a fit, so every one had a tolerance,
-- and the tolerances fought each other.
--
-- `LocalGrid` already solved it. Each block part is sampled on ITS OWN axes, so
-- a rotated part's rim falls on whole lattice lines of its own grid. There is no
-- staircase to fit away, and therefore nothing to fit: the outline of the cell
-- mask IS the outline, exactly, at whatever angle the part is at.
--
-- Two consequences worth stating plainly, because they are what the fitting
-- pipeline could never have:
--
--   * HEIGHTS ARE EXACT BY CONSTRUCTION. A vertex is a point on the part's own
--     face plane, computed as origin + u*uc + v*vc. Nothing is sampled from a
--     neighbouring cell and nothing is interpolated, so a vertex cannot pick up
--     a height from the wrong side of a cliff.
--   * THE RIM IS THE PART'S RIM. The cell lattice covers floor(2*ext/step)
--     whole cells, so the mask stops up to one step short of the real face
--     edge. `LocalGrid` keeps `center`, `uExt` and `vExt` for exactly this
--     reason: an outline vertex sitting on the mask's outermost lattice line is
--     mapped to the part's true extent instead. The polygon lands on the part's
--     face, not on the lattice that approximates it.
--
-- Non-block parts (Unions, MeshParts, wedges) have no meaningful surface axes,
-- so `LocalGrid` gives them a world-aligned fallback grid. Those are traced by
-- the same code in world XZ, and they are the only place a staircased rim can
-- still appear. That is a known limit of the fallback, not of this module.

local Boundary = {}

export type Config = {
	-- Fallback grids only. A world-aligned grid can hold two surfaces at once
	-- (a mesh with a ledge), so its trace still needs to know what counts as a
	-- cliff. Block grids are a single plane and never consult this.
	stepTol: number?,
}

local DEFAULT = {
	stepTol = 2.2,
}

local function merged(cfg): any
	local c = {}
	for k, v in pairs(DEFAULT) do c[k] = v end
	if cfg then for k, v in pairs(cfg) do if v ~= nil then c[k] = v end end end
	return c
end

--------------------------------------------------------------------------
-- Contour of a cell mask, in whatever integer lattice the mask is indexed by.
--
-- Boundary EDGES are chained rather than boundary cells walked: the boundary of
-- a set of cells is a closed loop by construction, so there is no open end to
-- chase and no tolerance involved.
--
-- `cliff(a, b)` reports that two lattice-adjacent cells are not continuous
-- ground even though both are in the mask. On a block grid it is never true.
--------------------------------------------------------------------------

local DIRS = { { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } }

-- lattice corners of cell (u,v)'s edge in direction di, oriented so the loop
-- runs with the interior on its left
local function corners(u: number, v: number, di: number)
	if di == 1 then return { u + 1, v }, { u + 1, v + 1 }
	elseif di == 2 then return { u + 1, v + 1 }, { u, v + 1 }
	elseif di == 3 then return { u, v + 1 }, { u, v }
	else return { u, v }, { u + 1, v } end
end

local function traceMask(cells: {any}, index: {[string]: any}, cliff: ((any, any) -> boolean)?)
	local segs: {any} = {}
	local byStart: {[string]: {any}} = {}
	for _, cell in ipairs(cells) do
		for di, d in ipairs(DIRS) do
			local nb = index[(cell.ui + d[1]) .. ":" .. (cell.vi + d[2])]
			if not nb or (cliff and cliff(cell, nb)) then
				local a, b = corners(cell.ui, cell.vi, di)
				local s = { a = a, b = b, cell = cell }
				segs[#segs + 1] = s
				local k = a[1] .. "," .. a[2]
				local bucket = byStart[k]
				if not bucket then bucket = {}; byStart[k] = bucket end
				bucket[#bucket + 1] = s
			end
		end
	end

	-- CHAINING AT A JUNCTION is not a free choice. A cliff between two cells
	-- that are both in the mask emits an edge from EACH side, on the same
	-- lattice edge in opposite directions, because the upper rim and the lower
	-- rim really are two different boundaries that coincide in plan. The mask's
	-- outline therefore has zero-width slits whose ends are lattice vertices
	-- with four boundary edges leaving them. Taking whichever candidate came
	-- first there splices one rim onto the other and drops degenerate slivers.
	--
	-- The rule that makes it deterministic is the standard face traversal of an
	-- embedded planar graph: arriving along an edge, leave on the next edge
	-- CLOCKWISE from it. With the interior kept on the left by `corners`, that
	-- walks each rim whole and turns a slit around at its tip.
	local function angle(s: any): number
		return math.atan2(s.b[2] - s.a[2], s.b[1] - s.a[1])
	end
	local TAU = math.pi * 2
	local used: {[any]: boolean} = {}
	local loops: {{any}} = {}
	for _, s0 in ipairs(segs) do
		if not used[s0] then
			local loop, cur = {}, s0
			while cur and not used[cur] do
				used[cur] = true
				loop[#loop + 1] = cur
				local cands = byStart[cur.b[1] .. "," .. cur.b[2]]
				local nxt, bestTurn = nil, math.huge
				if cands then
					local back = angle(cur) + math.pi
					for _, s in ipairs(cands) do
						if not used[s] then
							local turn = (back - angle(s)) % TAU
							-- straight back is the full turn, so it is taken only
							-- when nothing else is left: the tip of a slit
							if turn <= 1e-9 then turn = TAU end
							if turn < bestTurn then nxt, bestTurn = s, turn end
						end
					end
				end
				cur = nxt
			end
			if #loop >= 4 then loops[#loops + 1] = loop end
		end
	end
	return loops
end

-- The ordered lattice corners of a loop, with the duplicate a corner cell
-- contributes collapsed, and with runs along one lattice line reduced to their
-- endpoints. On a local grid a straight rim is exactly one lattice line, so
-- this is not a fit and has no tolerance -- it is just not repeating a point.
local function loopCorners(loop: {any}): {{number}}
	local pts: {{number}} = {}
	for _, s in ipairs(loop) do
		local last = pts[#pts]
		if not last or last[1] ~= s.a[1] or last[2] ~= s.a[2] then
			pts[#pts + 1] = { s.a[1], s.a[2] }
		end
	end
	if #pts > 1 and pts[1][1] == pts[#pts][1] and pts[1][2] == pts[#pts][2] then
		pts[#pts] = nil
	end
	if #pts < 3 then return {} end
	local out: {{number}} = {}
	for i = 1, #pts do
		local p = pts[i]
		local a = pts[(i - 2) % #pts + 1]
		local b = pts[i % #pts + 1]
		-- keep p only where the direction changes
		local d1u, d1v = p[1] - a[1], p[2] - a[2]
		local d2u, d2v = b[1] - p[1], b[2] - p[2]
		if d1u * d2v - d1v * d2u ~= 0 then out[#out + 1] = p end
	end
	return #out >= 3 and out or pts
end

local function signedArea(pts: {Vector3}): number
	local a2 = 0
	for i = 1, #pts do
		local p, q = pts[i], pts[i % #pts + 1]
		a2 += p.X * q.Z - q.X * p.Z
	end
	return a2 * 0.5
end

--------------------------------------------------------------------------
-- One grid -> its rings, in world space
--------------------------------------------------------------------------

local function ringsOfGrid(g: any, c: any)
	local isBlock = not g.fallback and g.n ~= nil
	local step = g.step

	local toWorld
	if isBlock then
		-- Index -> local face coordinate, RIM-EXACT at the extremes. The lattice
		-- covers nu whole cells of `step`, which is up to one step short of the
		-- face's real half-extent, so the outermost lattice line is mapped to the
		-- part's own extent instead of to the lattice. This is the whole reason
		-- LocalGrid keeps center/uExt/vExt.
		local nu = math.max(1, math.floor(2 * g.uExt / step + 1e-6))
		local nv = math.max(1, math.floor(2 * g.vExt / step + 1e-6))
		toWorld = function(iu: number, iv: number): Vector3
			local uc = (iu <= 0) and 0 or (iu >= nu and 2 * g.uExt or iu * step)
			local vc = (iv <= 0) and 0 or (iv >= nv and 2 * g.vExt or iv * step)
			return g.origin + g.u * uc + g.v * vc
		end
	else
		-- Fallback grid: ui/vi are world XZ lattice indices, and the surface is
		-- not one plane, so a corner takes the height of the nearest cell that
		-- owns it. No search and no fit -- the four cells touching a lattice
		-- corner are the only candidates there are.
		toWorld = function(iu: number, iv: number): Vector3
			local best, bestY = nil, 0
			for du = -1, 0 do
				for dv = -1, 0 do
					local cell = g.index[(iu + du) .. ":" .. (iv + dv)]
					if cell and (not best or cell.pos.Y > bestY) then best, bestY = cell, cell.pos.Y end
				end
			end
			return Vector3.new(iu * step, best and best.pos.Y or 0, iv * step)
		end
	end

	local cliff = nil
	if not isBlock then
		cliff = function(a: any, b: any): boolean
			return math.abs(a.pos.Y - b.pos.Y) > c.stepTol
		end
	end

	local rings: {any} = {}
	for _, loop in ipairs(traceMask(g.cells, g.index, cliff)) do
		local lat = loopCorners(loop)
		if #lat >= 3 then
			local verts = table.create(#lat)
			for i, p in ipairs(lat) do verts[i] = toWorld(p[1], p[2]) end
			rings[#rings + 1] = { verts = verts, area = signedArea(verts) }
		end
	end
	return rings
end

--------------------------------------------------------------------------
-- Boundary.fromLocal — the entry point
--------------------------------------------------------------------------

function Boundary.fromLocal(localData: any, cfg: Config?)
	local c = merged(cfg)
	c.stepTol = c.stepTol or 2.2
	local t0 = os.clock()

	local regions: {any} = {}
	local stats = {
		parts = 0, block = 0, fallback = 0,
		regions = 0, holes = 0, verts = 0, emptyGrids = 0,
	}

	for part, g in pairs(localData.grids) do
		stats.parts += 1
		if g.fallback then stats.fallback += 1 else stats.block += 1 end
		local rings = ringsOfGrid(g, c)
		if #rings == 0 then
			stats.emptyGrids += 1
			continue
		end
		-- Largest magnitude is the outer ring. Magnitude, never the sign: this
		-- project has been bitten by a handedness assumption before.
		local outer = rings[1]
		for _, r in ipairs(rings) do
			if math.abs(r.area) > math.abs(outer.area) then outer = r end
		end
		local region = {
			part = part,
			fallback = g.fallback,
			verts = outer.verts,
			area = math.abs(outer.area),
			cells = #g.cells,
			dead = #g.dead,
			holes = {},
		}
		for _, r in ipairs(rings) do
			if r ~= outer then
				region.holes[#region.holes + 1] = { verts = r.verts, area = math.abs(r.area) }
			end
		end
		stats.verts += #region.verts
		for _, h in ipairs(region.holes) do stats.verts += #h.verts end
		stats.regions += 1
		stats.holes += #region.holes
		regions[#regions + 1] = region
	end

	stats.seconds = os.clock() - t0
	return { regions = regions, stats = stats, config = c }
end

--------------------------------------------------------------------------

function Boundary.visualize(res: any, parent: Instance?)
	local root = parent or workspace
	local old = root:FindFirstChild("NavGen_Boundary")
	if old then old:Destroy() end
	local folder = Instance.new("Folder")
	folder.Name = "NavGen_Boundary"
	folder.Parent = root

	local function bar(a: Vector3, b: Vector3, colour: Color3, thick: number, into: Instance)
		local d = b - a
		if d.Magnitude < 1e-4 then return end
		local p = Instance.new("Part")
		p.Anchored = true; p.CanCollide = false; p.CanQuery = false; p.CanTouch = false
		p.Material = Enum.Material.Neon
		p.Color = colour
		p.Size = Vector3.new(thick, thick, d.Magnitude)
		p.CFrame = CFrame.lookAt(a + d * 0.5, b) + Vector3.new(0, 0.08, 0)
		p.Parent = into
	end

	-- Outer rings red, holes cyan, and a FALLBACK part's outline in amber --
	-- those are the world-aligned ones, the only place a staircased rim can
	-- still appear, so they should be identifiable without counting anything.
	for ri, r in ipairs(res.regions) do
		local sub = Instance.new("Folder")
		sub.Name = string.format("R%d_%s_c%d_h%d%s", ri, r.part.Name, r.cells, #r.holes,
			r.fallback and "_FALLBACK" or "")
		sub.Parent = folder
		local outerCol = r.fallback and Color3.fromRGB(255, 170, 60) or Color3.fromRGB(255, 80, 80)
		local v = r.verts
		for i = 1, #v do
			bar(v[i], v[i % #v + 1], outerCol, 0.18, sub)
		end
		for _, hr in ipairs(r.holes) do
			local hv = hr.verts
			for i = 1, #hv do
				bar(hv[i], hv[i % #hv + 1], Color3.fromRGB(80, 200, 255), 0.18, sub)
			end
		end
	end
	return folder
end

return Boundary
